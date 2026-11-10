import pandas as pd, numpy as np
PB="/well/blagden/users/eic789/PerturbBridge"
H=pd.read_parquet(f"{PB}/14_inverse/dict_splicing_rbp_HepG2.parquet"); K=pd.read_parquet(f"{PB}/14_inverse/dict_splicing_rbp_K562.parquet")
print("shapes", H.shape, K.shape, flush=True)
H=H.pivot(index="target",columns="event_id",values="dpsi_cc"); K=K.pivot(index="target",columns="event_id",values="dpsi_cc")
ev=H.columns.intersection(K.columns); rb=H.index.intersection(K.index)
ok=(H.loc[rb,ev].notna().mean(0).values>0.8)&(K.loc[rb,ev].notna().mean(0).values>0.8); ev=ev[ok]
A=H.loc[rb,ev].fillna(0).values.astype(np.float32); B=K.loc[rb,ev].fillna(0).values.astype(np.float32)
print("shared RBPs",len(rb),"shared well-covered events",len(ev), flush=True)
def rowcorr(X,Y):
    Xc=X-X.mean(1,keepdims=True); Yc=Y-Y.mean(1,keepdims=True); return (Xc*Yc).sum(1)/np.sqrt((Xc**2).sum(1)*(Yc**2).sum(1)+1e-12)
rc=rowcorr(A,B); print("event-level cross-cell r per RBP: median %.3f, q75 %.3f, q90 %.3f"%(np.median(rc),np.quantile(rc,.75),np.quantile(rc,.9)), flush=True)
rng=np.random.default_rng(0); R=[]; V=[]
for rep in range(10):
    p=rng.permutation(len(rb)); tr,te=p[:len(rb)//2],p[len(rb)//2:]
    mu=B[tr].mean(0); U,S,Vt=np.linalg.svd(B[tr]-mu,full_matrices=False); k=10
    sK=(B[te]-mu)@Vt[:k].T; sH=(A[te]-A[tr].mean(0))@Vt[:k].T
    R.append([np.corrcoef(sK[:,i],sH[:,i])[0,1] for i in range(k)]); V.append(S[:k]**2/(S**2).sum())
R=np.array(R).mean(0); V=np.array(V).mean(0)
print("factors fit on K562 with half of the RBPs; scores of the other half compared across cell lines")
for i in range(10): print(" PC%d  var %.3f  cross-cell r of scores %.2f"%(i+1,V[i],R[i]))
print("cum var top10 %.3f"%V.sum())
