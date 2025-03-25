import numpy as np
import cmath
import math
from scipy.integrate import quad

def heston_char_func(z, S0, r, T, v0, theta, gamma, sigma, rho):
    """Calculate the characteristic function of the Heston model"""
    i = 1j  # complex unit
    xi = gamma - rho * sigma * i * z
    d = cmath.sqrt((rho * sigma * i * z - gamma)**2 + sigma**2 * z * (i + z))
    g = (xi - d) / (xi + d)
    
    # Avoid numerical issues caused by g
    if abs(1 - g) < 1e-10:
        g += 1e-10
    
    C = (gamma * theta / sigma**2) * (
        (xi - d) * T - 2 * cmath.log((1 - g * cmath.exp(-d * T)) / (1 - g))
    )
    D = ((xi - d) / sigma**2) * ((1 - cmath.exp(-d * T)) / (1 - g * cmath.exp(-d * T)))
    phi = cmath.exp(i * z * cmath.log(S0) + i * z * r * T + C + D * v0)
    return phi

def heston_call_price(S0, K, T, r, v0, theta, gamma, sigma, rho, upper_limit=1000):
    """Calculate European call option price using the Heston model"""
    # Define integration functions
    def f1(u):
        phi_u_minus_i = heston_char_func(u - 1j, S0, r, T, v0, theta, gamma, sigma, rho)
        integrand = cmath.exp(-1j * u * cmath.log(K)) * phi_u_minus_i / (1j * u * S0 * math.exp(r * T))
        return integrand.real

    def f2(u):
        phi_u = heston_char_func(u, S0, r, T, v0, theta, gamma, sigma, rho)
        integrand = cmath.exp(-1j * u * cmath.log(K)) * phi_u / (1j * u)
        return integrand.real

    # Numerical integration, avoiding singularity at u=0
    P1_integral, _ = quad(f1, 1e-6, upper_limit)
    P2_integral, _ = quad(f2, 1e-6, upper_limit)

    P1 = 0.5 + (1 / math.pi) * P1_integral
    P2 = 0.5 + (1 / math.pi) * P2_integral

    # Calculate call option price
    call_price = S0 * P1 - K * math.exp(-r * T) * P2
    return call_price

def heston_put_price(S0, K, T, r, v0, theta, gamma, sigma, rho, upper_limit=1000):
    """Calculate European put option price using the Heston model"""
    # Call the call option price function
    call_price = heston_call_price(S0, K, T, r, v0, theta, gamma, sigma, rho, upper_limit)
    # Use put-call parity to calculate put option price
    put_price = call_price - S0 + K * math.exp(-r * T)
    return put_price

# Test parameters
S0 = 100    # Current stock price
strike_prices = [80, 90, 100, 110, 120]  # Multiple strike prices
T = 1       # Time to maturity (years)
r = 0    # Risk-free rate
theta = 0.3156 # Long-term variance
v0 = 0.0392   # Initial variance
gamma = 0.1  # Mean reversion speed
nu = 0.331
sigma = gamma*nu  # Volatility of volatility
rho = -0.681   # Correlation coefficient

# Calculate call and put option prices for different strike prices
call_prices = [heston_call_price(S0, K, T, r, v0, theta, gamma, sigma, rho) for K in strike_prices]
put_prices = [heston_put_price(S0, K, T, r, v0, theta, gamma, sigma, rho) for K in strike_prices]

# Print results table
print("Heston Model European Call and Put Option Prices:")
print("-" * 50)
print(f"{'Strike':^10} | {'Call Price':^15} | {'Put Price':^15}")
print("-" * 50)
for K, call, put in zip(strike_prices, call_prices, put_prices):
    print(f"{K:^10} | {call:^15.4f} | {put:^15.4f}")
print("-" * 50)