from math import ceil, log, floor

p = 757
q = 7879
w = 242
I = 128
tau = 4
tau0 = 3011
tau1 = 33
tau2 = 1995
tau3 = 1978

# From reference sage script
def Encode(R,M):
  if len(M) == 0: return []
  S = []
  if len(M) == 1:
    r,m = R[0],M[0]
    while m > 1:
      S += [r%256]
      r,m = r//256,(m+255)//256
    return S
  R2,M2 = [],[]
  for i in range(0,len(M)-1,2):
    m,r = M[i]*M[i+1],R[i]+M[i]*R[i+1]
    while m >= 16384:
      S += [r%256]
      r,m = r//256,(m+255)//256
    R2 += [r]
    M2 += [m]
  if len(M)&1:
    R2 += [R[-1]]
    M2 += [M[-1]]
  return S+Encode(R2,M2)

def barrett(p, q, w, inval=0):
    if inval == 0:
        inval = 3*w*(q-1)//2
    print()
    print(f' p={p} q={q} w={w}')
    print(f"Unreduced coefficients ≤2^{log(inval,2):.2f} = {inval/q:.1f}q")
    best = (0,0,0)
    for s in range(2,32):
        num = int(floor(2**s/q+.5))
        if not num:
            continue
        err = abs(1/q - num/(2**s))
        max_in = min(2**31 / num, 1/(q*err))
        if max_in > best[0]:
            best = (max_in, num, s)
    max_in, num, s = best
    print(f'Barrett reduction {num}x/2^{s} requires |x|≤{max_in/q:.2f}q')
    print("Initial reductions:")
    ok = False
    for s in range(2,32):
        num = int(floor(2**s/q+.5))
        if not num:
            continue
        err = abs(1/q - num/(2**s))
        res = ceil(inval * err - 1/q)+1
        if res*q > max_in:
            continue
        if inval > 2**31 / num:
            continue
        print(f'  |{num}x/2^{s}| ≤ {res}q')
        ok = True
    if ok:
        return True
    print( '  <none found>')
    return False

def full_barrett(q):
    print()
    print(f' Full Barrett reduction for q={q}')
    best = (0,0,0)
    for s in range(2,63):
        num = int(floor(2**s/q+.5))
        if not num:
            continue
        err = abs(1/q - num/(2**s))
        max_in = min(2**64 / num, 1/(q*err))
        if max_in > best[0]:
            best = (max_in, num, s)
    max_in, num, s = best
    print(f'Barrett reduction {num}x/2^{s} requires |x|≤2^{log(max_in,2):.1f}')


print("Size of rounded elements: ", len(Encode([0]*p, [(q-1)//3+1]*p)))
print("Size of short element:", int(ceil(p/4)))
print("Size of ct1:", int(ceil((I * int(log(tau,2)))/8)))

barrett(p, q, w)

full_barrett(q)

print("Checking rounding …")
for x in range(-2**14, 2**14):
    assert abs(3* ((10923*x + 16384) >> 15)-x) <= 1


pp = 2**9*3
qq = 956929
zeta = 2411

assert qq > w*q//2

for i in range(1, pp):
    assert pow(zeta, i, qq) != 1
assert pow(zeta, pp, qq) == 1

print(pow(zeta, pp//2, qq))
print(qq)


import gmpy2
q3 = w*(q-1)//2
print(f'minval= {q3}')
print(f'minval2= {w*(q-1)//2}, {log(w*(q-1)//2,2)}')
print(f'q2^9= {q*2**9}')
print(f'2*p={2*p}')
print(q)
while True:
    q3 = gmpy2.next_prime(q3)
    if (q3 - 1) % 2**9 == 0:
        print(q3)
        break
