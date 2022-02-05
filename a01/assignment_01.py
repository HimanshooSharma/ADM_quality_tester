import numpy as np
from numpy.linalg import matrix_power


def calculate_nth_step(transition_matrix, pi_zero, n):
    """A function that calculates probability vector after step n

    Args:
        transition_matrix (2D numpy array): The state transition matrix
        pi_zero (1D numpy array): The initial probability vector
        n (int): The number of steps

    Returns:
        1D numpy array: Probability vector at step n
    """

    return matrix_power(transition_matrix, n).dot(pi_zero)


def find_steady_state(transition_matrix, pi):
    """A function that finds the steady state of the DTMC

    Args:
        transition_matrix (2D numpy array): The state transition matrix
        pi_zero (1D numpy array): The initial probability vector]

    Returns:
        1D numpy array, int: The steady state probability vector, no. of steps to reach the steady state
    """

    next_pi = None
    step = 0

    while True:
        next_pi = transition_matrix.dot(pi)
        if (pi == next_pi).all():
            break
        pi = next_pi
        step += 1

    return next_pi, step


transition_probability_matrix = np.array([[0.6, 0.4], [0.3, 0.7]])

initial_pi = np.array([1, 0])  ## Source 0 is active

pi_after_8_steps = calculate_nth_step(transition_probability_matrix, initial_pi, 8)

print(f"Probability that source 0 is active after 8 minutes is {pi_after_8_steps[0]}")

pi_after_9_steps = calculate_nth_step(transition_probability_matrix, initial_pi, 9)

prob_producing_ok_at_9_steps = (pi_after_9_steps[0] * 0.9) + (
    pi_after_9_steps[1] * 0.95
)


print(
    f"Probability of producing an OK item in the next minute: {prob_producing_ok_at_9_steps}"
)

steady_pi, steady_state_steps = find_steady_state(
    transition_probability_matrix, initial_pi
)

print(f"Stationary solution reached at {steady_state_steps} steps")

print(f"Probability for source 0 to be in steady state : {steady_pi[0]}")

prob_OK_in_steady = (steady_pi[0] * 0.9) + (steady_pi[1] * 0.95)

print(
    f"Average probability of producing an OK item in steady state : {prob_OK_in_steady}"
)
