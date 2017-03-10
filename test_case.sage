from ntru_attack import attack

for i in range(30, 80, 10):
    for _ in range(4):
        attack(m=256, q = next_prime(2^i), r = 4, sigma=2, subfield_only=True)

for i in range(30, 80, 10):
    for _ in range(4):
        attack(m=256, q = next_prime(2^i), r = 2, sigma=2, subfield_only=True)

