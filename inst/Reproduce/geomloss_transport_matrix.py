import torch
import pykeops
import geomloss
import numpy as np

torch.manual_seed(23432)

def round_pi(raw_pi, a, b):
  n = len(a)
  m = len(b)
  x = a / torch.sum(raw_pi, 1)
  for i in range(n):
    if(x[i] > 1):
      x[i] = 1
  
  y = b / torch.sum(raw_pi, 0)
  for i in range(m):
    if(y[i] > 1):
      y[i] = 1
  
  pi_prime = (x.reshape(-1,1) * raw_pi) * y
  err_row = a - pi_prime.sum(1)
  err_col = b - pi_prime.sum(0)
  
  return(pi_prime + (err_row.reshape(-1,1) * err_col / err_row.abs().sum()) )


# initalize data
x = torch.randn(100,5)
y = torch.randn(1000,5)
a = torch.FloatTensor([1/100] * 100)
b = torch.FloatTensor([1/1000] * 1000)

blur_list = [0.05,.5,5,50,500] # lambda values
dual_loss = [None] * len(blur_list)
primal_loss = [None] * len(blur_list)
gamma_log = torch.zeros(100,1000)

for i in range(len(blur_list)):
  blur = blur_list[i] # set blur
  Sinkhorn = geomloss.SamplesLoss("sinkhorn", p = 2, blur = blur, 
    cost = geomloss.utils.squared_distances,
                      debias = False,
                      potentials = True,
                      verbose = True,
                      scaling = 0.99)
  output = Sinkhorn(a, x, b, y)
  
  # get potentials
  f = output[0]
  g = output[1]
  
  #calculate loss using implied transportation matrix
  cost = torch.cdist(x,y, p = 2)**2
  gamma = round_pi(((f.reshape(-1,1) + g - cost)/blur).exp() * a.reshape(-1,1) * b, a,b)
  # gamma = ((f.reshape(-1,1) + g - cost)/blur).exp() * a.reshape(-1,1) * b
  for j in range(1000):
    for k in range(100):
      if(gamma[k,j] == 0.0):
        gamma_log[k,j] = 0.0
      else :
        gamma_log[k,j] = (gamma[k,j]).log()
  primal_loss[i] = torch.sum(gamma * cost) - blur * (gamma *( gamma_log - (a.reshape(-1,1)*b).log())).sum()
  
  # compare to dual loss
  dual_loss[i] = torch.sum(f * a) + torch.sum(g * b)

print(primal_loss)
print(dual_loss)
# values don't align
# [tensor(0.0247), tensor(0.3793), tensor(10.3293), tensor(10.0444), tensor(10.3681)]
# tensor(1.7243), tensor(2.7584), tensor(9.9867), tensor(10.4062), tensor(10.3986)]
