import numpy as np
from numpy.linalg import matrix_power


def calculate_n_steps_with_discretization(Q, pi_zero, n, delta_t):

    # Intializing the identity matrix I
    I = np.array([[1, 0], [0, 1]])

    # calculating transition matrix as delta_t*Q + I
    transition_matrix = np.multiply(delta_t, Q) + I

    return matrix_power(transition_matrix, n).dot(pi_zero)


def calculate_steady_state(Q, pi, delta_t):

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


generator_matrix = np.array([[-0.33, 0.5], [0.33, -0.5]])

initial_pi = np.array([1, 0])  ## Source 0 is active

print(
    f"Probability that source 0 is active after 8 minutes for delta_t 2 : {calculate_n_steps_with_discretization(generator_matrix,initial_pi,8,2)[0]}"
)
print(
    f"Probability that source 0 is active after 8 minutes for delta_t 1 : {calculate_n_steps_with_discretization(generator_matrix,initial_pi,8,1)[0]}"
)
print(
    f"Probability that source 0 is active after 8 minutes for delta_t 0.5 : {calculate_n_steps_with_discretization(generator_matrix,initial_pi,8,0.5)[0]}"
)
print(
    f"Probability that source 0 is active after 8 minutes for delta_t 0.25 : {calculate_n_steps_with_discretization(generator_matrix,initial_pi,8,0.25)[0]}"
)
print(
    f"Probability that source 0 is active after 8 minutes for delta_t 0.1 : {calculate_n_steps_with_discretization(generator_matrix,initial_pi,8,0.1)[0]}"
)

print(
    f"Cannot find steady state for detla_t 2 with precision of over 15 decimal places"
)
print(
    f"Probability of source 0 being active in steady state for delta_t 1 : {calculate_steady_state(generator_matrix,initial_pi,1)[0]}"
)
print(
    f"Probability of source 0 being active in steady state for delta_t 0.5 : {calculate_steady_state(generator_matrix,initial_pi,0.5)[0]}"
)
print(
    f"Probability of source 0 being active in steady state for delta_t 0.25 : {calculate_steady_state(generator_matrix,initial_pi,0.25)[0]}"
)
print(
    f"Probability of source 0 being active in steady state for delta_t 0.1 : {calculate_steady_state(generator_matrix,initial_pi,0.1)[0]}"
)

steady_state_probability = calculate_steady_state(generator_matrix, initial_pi, 0.1)

average_OK_in_ss = (steady_state_probability[0] * 0.90) + (
    steady_state_probability[1] * 0.95
)

print(f"Average probability of testing an item OK in steady state : {average_OK_in_ss}")
