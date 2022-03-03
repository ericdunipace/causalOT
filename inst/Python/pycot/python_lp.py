import torch

def lp_loss(x, y, p):
  D = torch.cdist(x, y, p=p, compute_mode='use_mm_for_euclid_dist_if_necessary')
  
  return D**p


def l1_loss(x, y):
  D = torch.cdist(x, y, p=1.0, compute_mode='use_mm_for_euclid_dist_if_necessary')
  return D
