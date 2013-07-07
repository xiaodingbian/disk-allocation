import logging
from numpy import *
#from pylab import *

class matrixUtil:
  def __init__(self, lines, columns):
    self.lines = lines
    self.columns = columns
    self.data = self.generate()

  def generate(self):
    randMatr = zeros((self.lines, self.columns))
    if self.lines > self.columns:
      for i in range(self.columns):
        randMatr[i][i] = 1
      for i in range(self.columns, self.lines):
        k = random.randint(0, self.columns)
        randMatr[i][k] = 1
    else:
      for i in range(self.lines):
        randMatr[i][i] = 1
      for i in range(self.lines, self.columns):
        k = random.randint(0, self.lines)
        randMatr[k][i] = 1

    for j in range(self.columns):
      weight = random.randint(1, int(self.lines * 0.5))
      for i in range(self.lines):
        if random.randint(1, self.lines) <= weight:
          randMatr[i][j] = 1
    return randMatr

  def steadyDistribute(self, loop, epison):
    distr = self.data.copy()
    for k in range(loop):
      self.normalize(distr, self.lines, self.columns)
      logging.debug("steady-state distribution at %d Loop", k)
      logging.debug(distr)
      done = self.validateAndAdjust(distr, self.lines, self.columns, 1.0 * self.lines / self.columns, epison)
      if done:
        logging.debug("Already accurate enough")
        break
      else:
        logging.debug("After adjust at %d Loop:", k)
        logging.debug(distr)

    self.normalize(distr, self.lines, self.columns)
    logging.info("final steady-state distribution:")
    logging.info(distr)
    self.A = distr

  def reArrange(self):
    self.B = []
    self.mapBtoA = []
    for i in range(self.lines):
      tempArray = []
      tempMap = {}
      for j in range(self.columns):
        if not (self.A[i][j] in tempMap):
          tempMap[self.A[i][j]] = []
        tempMap[self.A[i][j]].append(j)
      logging.debug("mapping for VM: %d", i)
      logging.debug(tempMap)
      for j in range(self.columns):
        if self.A[i][j] > 0:
          tempArray.append(self.A[i][j])
      tempArray.sort()
      b = zeros(len(tempArray))
      k = 0
      j = 0
      while j < len(tempArray):
        b[k] = tempArray[j]
        k += 1
        j += 2
      k = len(tempArray) - 1
      j = 1
      while j < len(tempArray):
        b[k] = tempArray[j]
        k -= 1
        j += 2
      logging.debug("Rearrange result for %d sub Array:")
      logging.debug(b)
      self.B.append(b)

      mapbtoa = zeros(len(b))
      for k in range(len(b)):
        mapbtoa[k] = tempMap[b[k]].pop()
      self.mapBtoA.append(mapbtoa)
    logging.info("After rearranged:")
    logging.info(self.B)
    logging.info("Mapping matrix:")
    logging.info(self.mapBtoA)


  def calTranMatrix(self):
    self.transB = []
    for t in range(self.lines):
      b = self.B[t]
      beta = zeros(len(b))
      for j in range(len(b)):
        beta[j] = 1
        start = (j - 1) % len(b)
        temp = b[start]
        while temp < b[j]:
          beta[j] += 1
          start = (start - 1) % len(b)
          temp += b[start]
      logging.debug("Beta for VM %d", t)
      logging.debug(beta)
      self.calTranMatrixForVM(t, beta, 100, 0.001)
  
  def calTranMatrixForVM(self, vmid, beta, loop, epison):
    b = self.B[vmid]
    size = len(b)
    Pb = self.initPb(beta, size)

    logging.info("initialize Pb for VM: %d", vmid)
    logging.info(Pb)

    ret = self.adjustTransient(vmid, Pb, size, loop, epison, beta)
    #logging.info("ret value: %d", ret)
    while ret != -1:
      logging.debug("old beta")
      logging.debug(beta)
      k = beta[0]
      index = 0
      for t in range(size):
        if k > beta[t]:
          k = beta[t]
          index = t
      if k == size:
        break
      beta[index] += 1

#      if beta[ret] < size:
#        beta[ret] += 1
#      else:
#        k = beta[0];
#        index = 0;
#        for t in range(size):
#          if k > beta[t]:
#            k = beta[t]
#            index = t
#        if k == size:
#          break
#        beta[index] += 1

      logging.debug("new beta")
      logging.debug(beta)
      newPb = self.initPb(beta, size)
      self.normalize(newPb, size, size)
      for i in range(size):
        for j in range(size):
          Pb[i][j] = (Pb[i][j] + newPb[i][j]) * 0.5
      ret = self.adjustTransient(vmid, Pb, size, loop, epison, beta)
    self.transB.append(Pb)

  def initPb(self, beta, size):
    Pb = zeros((size, size))
    for j in range(size):
      nozeros = int(beta[j])
      for i in range(nozeros):
        index = (j - i - 1) % size
        Pb[index][j] = 1
    return Pb


  def adjustTransient(self, vmid, trans, size, loop, epison, beta):

    tempb = zeros(size)
    for k in range(loop):
      self.normalize(trans, size, size)
      b = self.B[vmid]

      for i in range(size):
        tempb[i] = 0
        for j in range(size):
          tempb[i] += b[j] * trans[j][i]

      error = 0
      for i in range(size):
        error += (b[i] - tempb[i]) * (b[i] - tempb[i])

      logging.debug("the error: %f", error)
      if error < epison:
        logging.info("Final Transient matrix:")
        logging.info(trans)
        return -1

      for i in range(size):
        for j in range(size):
          trans[i][j] = trans[i][j] * b[j] / tempb[j]

    self.normalize(trans, size, size)
    index = 0
    value = 0
    for i in range(size):
      if tempb[i] / b[i] > 1 and ( index == 0 or beta[i] < value):
        value = beta[i]
        index = i
    logging.debug("Best Transient matrix under current template:")
    logging.debug(trans)
    return index

  def normalize(self, array, n, m):
    for i in range(n):
      sum = 0
      for j in range(m):
        sum += array[i][j]
      factor = 1.0 / sum
      for j in range(m):
        array[i][j] = array[i][j] * factor

  def validateAndAdjust(self, array, n, m, refValue, epison):
    ret = True
    sum = zeros(m)
    for j in range(m):
      sum[j] = 0
      for i in range(n):
        sum[j] += array[i][j]
      if abs(sum[j] - refValue) > epison:
        ret = False
    logging.debug("Total of each column:")
    logging.debug(sum)

    if ret:
      return ret

    for j in range(m):
      for i in range(n):
        array[i][j] = array[i][j] * 1.0 * n / m / sum[j]
    return ret

  def generateSeq(self):
    self.sequence = []
    for i in range(self.lines):
      b = self.B[i]
      trans = self.transB[i]
      logging.info("b matrix for VM %d:", i)
      logging.info(b)
      logging.info("Its transient matrix:")
      logging.info(trans)
      p = zeros(len(b))
      temp = random.randint(0, len(b))
      p[temp] = 1
      seq = []
      for j in range(1000):
        disk = self.selectDisk(p, trans, b)
        seq.append(disk)
      logging.info("sequence for VM %d:", i)
      logging.info(seq)
      self.sequence.append(seq)

  def evalByCos(self):
    for i in range(self.lines):
      logging.info("cos evaluation for VM %d:", i)
      T = 0
      S1 = 0
      S2 = 0
      index = 0
      b = self.B[i]
      fragment = len(b) * 50
      seq = self.sequence[i]
      actual = zeros(len(b))

      cosArray = []
      while index + fragment < len(seq):
        for j in range(fragment):
          actual[seq[index+j]] += 1.0 / fragment
        index += fragment
        total = 0
        for k in range(len(b)):
          T += actual[k] * b[k]
          S1 += actual[k] * actual[k]
          S2 += b[k] * b[k]
        logging.info(abs(T) / math.sqrt(S1 * S2))
        cosArray.append(abs(T) / math.sqrt(S1 * S2))
#      x = arange(0, len(cosArray), 1)
#      y = cosArray
#      plot(x, y, linewidth=1.0)
#    xlabel('the x')
#    ylabel('the cos value')
#    title('evaluation')
#    grid(True)
#    savefig("cos.png")
#    show()
          
  def selectDisk(self, p, trans, b):
    rand = random.uniform(0,1)
    index = 0
    sum = 0
    for index in range(len(b)):
      if index == len(b) - 1:
        break
      sum += b[index]
      if sum > rand:
        break
      else:
        index += 1
    temp = zeros(len(p))
    for i in range(len(p)):
      temp[i] = 0
      for j in range(len(p)):
        temp[i] += p[j] * trans[j][i]
    p = temp
    return index


