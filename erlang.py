import math
def erlangc(lam, mu, c):
    a = lam / mu; rho = a / c
    s = sum(a**k / math.factorial(k) for k in range(c)); top = a**c / (math.factorial(c) * (1 - rho))
    return top / (s + top)
def mmc(lam, mu, c, T):
    pw = erlangc(lam, mu, c); al = c * mu - lam
    wq = pw / al; R = wq + 1 / mu
    surv_wait = math.exp(-mu*T)*(1+mu*T) if abs(al-mu) < 1e-12 else (mu*math.exp(-al*T) - al*math.exp(-mu*T)) / (mu - al)
    p_gt = (1 - pw) * math.exp(-mu*T) + pw * surv_wait
    return dict(pw=pw, wq_ms=1000*wq, R_ms=1000*R, attain=1 - p_gt)
if __name__ == "__main__":
    mu = 15
    for lam in (40, 60, 80, 90, 110):
        r = mmc(lam/2, mu, 4, 0.25); p = mmc(lam, mu, 8, 0.25)
        print(lam, round(lam/120, 3), round(r['R_ms'], 1), round(p['R_ms'], 1), round(100*r['attain'], 1), round(100*p['attain'], 1), "wq", round(r['wq_ms'],1), round(p['wq_ms'],1))
    p = mmc(90, mu, 8, 0.22); print("pooled +30ms: mean", round(mmc(90,mu,8,.25)['R_ms']+30,1), "attain@220", round(100*p['attain'],1))
    print("break-even", round(mmc(45,mu,4,.25)['R_ms']-mmc(90,mu,8,.25)['R_ms'],1))
    print("170ms budget pooled", round(100*mmc(90,mu,8,0.17)['attain'],1))
    print("goodput", round(90*mmc(45,mu,4,.25)['attain'],1), round(90*mmc(90,mu,8,.25)['attain'],1))
    for lam in (60,90,120,150): h = math.ceil(lam/(0.75*mu)); print(lam, h, 2*h, round(lam/(h*mu),3))
    print("72 rps:", math.ceil(72/11.25))
