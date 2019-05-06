import torch
import math
def C60(K):
    lnZ = 0
    Tbond = torch.exp(torch.tensor([[-K,K],[K,-K]],dtype=torch.float64))
    TId = torch.tensor([[1,0],[0,1]],dtype=torch.float64)
    T = torch.einsum('ai,ib,ic->abc',(Tbond,Tbond,Tbond))
    T_max = torch.max(T)
    T = T/T_max
    lnZ += torch.log(T_max)*10
    T = torch.einsum('abc,ib,be,cj,cf,ef,eg,hf->aijgh',(T,TId,Tbond,Tbond,Tbond,Tbond,TId,Tbond))
    T_max = torch.max(T)
    T = T/T_max
    lnZ += torch.log(T_max)*10
    T = torch.einsum('aijgh,zmnhf->azijgfmn',(T,T))
    T_max = torch.max(T)
    T = T/T_max
    lnZ += torch.log(T_max)*5
    T = torch.einsum('azijgfmn,pa,aq,rz,zs->pigmrqjfns',(T,TId,Tbond,TId,Tbond))
    T_max = torch.max(T)
    T = T/T_max
    lnZ += torch.log(T_max)*5
    dim = T.shape[0]*T.shape[1]*T.shape[2]*T.shape[3]*T.shape[4]
    T = T.contiguous().view(dim,dim)
    T = torch.trace(T.mm(T).mm(T).mm(T).mm(T))
    lnZ += torch.log(T)
    return lnZ.item()/60

if __name__ == "__main__":
    C_T1 = C60(1.0)
    print('lnZ=',C_T1)
    C_T100 = 60*C60(100.0)
    C_T90 = 60*C60(90.0)
    E = -(C_T100-C_T90)/10
    N = math.exp(C_T100+100*E)
    print(N)
    