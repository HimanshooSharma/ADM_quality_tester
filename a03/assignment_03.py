import numpy as np
from numpy.linalg import matrix_power


def calculate_steady_state(Q, pi, delta_t):
    """A function to calculate the steady state of ctmc

    Args:
        Q (2D numpy array)): The genrator matrix
        pi (1D numpy array): The probability vector
        delta_t (int): The discreetization step

    Returns:
        1D numpy arrray: Probability vector in steady state
    """

    # Intializing the identity matrix I
    I = np.array([[1, 0], [0, 1]])

    # calculating transition matrix as delta_t*Q + I
    transition_matrix = np.multiply(delta_t, Q) + I

    next_pi = None

    while True:
        next_pi = transition_matrix.dot(pi)
        if (pi == next_pi).all():
            break
        pi = next_pi

    return next_pi


generator_matrix = np.array([[-0.1, 0.143], [0.1, -0.143]])

initial_pi = np.array([1, 0])  ## Source 0 is active as it has the token

print(
    f"Steady state vector for deta_t 2 : {calculate_steady_state(generator_matrix, initial_pi, 2)}"
)
print(
    f"Steady state vector for deta_t 1 : {calculate_steady_state(generator_matrix, initial_pi, 1)}"
)
print(
    f"Steady state vector for deta_t 0.5 : {calculate_steady_state(generator_matrix, initial_pi, 0.5)}"
)
print(
    f"Steady state vector for deta_t 0.25 : {calculate_steady_state(generator_matrix, initial_pi, 0.25)}"
)
print(
    f"Steady state vector for deta_t 0.1 : {calculate_steady_state(generator_matrix, initial_pi, 0.1)}"
)

steady_state_vector = calculate_steady_state(generator_matrix, initial_pi, 0.1)

troughput_s1_in_steady = steady_state_vector[1] * 0.95

print(f"Throughput of source 1 in steady state : {troughput_s1_in_steady}")
