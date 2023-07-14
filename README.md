# BistableFlipFlop
Code to run decision-making dynamics in a mutual-inhibition bistable attractor system

To use this package,

1. Install the [Julia language](https://julialang.org/).
2. Change directory to this package's directory
3. Open a Julia session
4. Within Julia, type `using Pkg; pkg"add PyPlot; add LinearAlgebra; add Statistics`
5. Within Julia, type `include("moduleBistableFlipFlop")`
6. You're set

From here, within the Julia session, you can type `decisionMake()` to run the main decision-making function.

# Assignment Questions

Here are the questions you should answer by running the `decisionMake()` function (defined below) different times with different parameter values (see below for a description of `decisionMake()` and how to run it with parameter values other than the defaults). But mostly I want you to have fun with this, and explore what happens as you change parameters, and think about how that might be interpreted as a model of the biology.

For each of these questions, think and comment on your answer.

1.    What happens to the reaction time and to the probability of `r1` being the winner as you change coherence `c`? Can you plot psychometric and reaction time functions?
2.    What happens as you change the overall excitation `mu0`? What if `mu0` is very small or negative? Can you interpret the results you see in biological terms? What if `mu0` is very very large? Can you interpret *those results* in biological terms?
3. What happens if you run with `noise_sigma=0, c=0`? What if you now make the noise really small, but non-zero, e.g., `noise_sigma=0.001`? What happens at very low coherence when the noise is very small, e.g., `noise_sigma=0.001, c=1`? Biologically, does it make sense that noise would need to be large in the real brain?
4.    Do you see a difference between error and correct reaction times? If "yes", can you speculate as to what causes that difference?
4.    For the default value of `wR=0`, no self-excitation, we can interpret each of the nullclines as the input-output function of the corresponding node, as described in lecture. For example, the blue line can be thought of as "given a certain level of activity of `r1`, this is the stead-state value of `r2`". Now run with much stronger self-excitation (and some stronger inhibition to balance that), as in `decisionMake(wR=6, wI=9)`. Is there still a one-to-one relationship between value of `r1` and value of `r2` on the blue line? What do you think happened? Think about a single neuron, feeding back excitation on itself (as in lecture), and what can happen if the positive feedback from the self-excitation is strong.
  
Again, most of all, have fun running it and exploring different parameters!

### Some references

Some relevant papers to read for this include

- Machens, Romo, and Brody, Machens, "Flexible Control of Mutual Inhibition," *Science* (2005)

- Wang, "Probabilistic decision making by slow reverberation in
cortical circuits," *Neuron* (2002)

- Wong and Wang, “A Recurrent Network Mechanism of Time Integration in Perceptual Decisions,” *J. Neurosci.* (2006)

*and even*

- Hopfield, “Neurons with Graded Response Have Collective Computational Properties like Those of Two-State Neurons,” *PNAS* (1984)


---



---



# The decisionMake() function

In the cell below, you'll run `decisionMake()`. This function runs bistable decision dynamics, following the mutual inhibition architecture described in class. Try running the cell right now to see what it produces-- `decisionMake()` can run as is, with nothing inside the parentheses, in which case all the parameters will take their default values. The function  runs a number of independent trials (by default 100 independent trials), all in parallel, each one with its own independent noise.

# Outputs of decisionMake()

**In the Main Axes**:

> Each trial starts from `r1=0, r2=0` in the main axes and leaves a little green trail as it goes; each big green dot is where one trial ended up by the end of the simulation.
>
> The **blue curve** is the `dr2/dt=0` nullcline, and the **red curve** is the `dr1/dt=0` nullcline. The **grey lines** are the decision boundaries; when a trajectory crosses the one on the right, an `r1=winner` decision is emitted, and when it crosses the upper one, an `r2=winner` decision is emitted.

The **smaller axes on left and bottom** are histograms of the reaction times (i.e., times when a decision was emitted.

When it is finished running, in addition to the plots it produces, `decisionMake()` prints out the value of

- `P1`: out of all the trials that reached a boundary, the fraction of them that reached the `r1=winner` boundary

- `ND`: the fraction of total trials run that had did not yet reached a decision boundary when the simulation finished.

- `RT1`: the average time to reach the bound for the `r1=winner` decisions

- `RT2`: the average time to reach the bound for the `r2=winner` decisions


# Changing the parameters of `decisionMake()`

`decisionMake()` uses a number of parameters that define how it runs. To change them from their default values, include the parameter you want to change within the parentheses, with an equal sign and the value you want it to take. For example, `c` stands for random dot coherence. To run at `c=40`, you would edit the cell below to now say

     decisionMake(c=40);

and then run it.  You can include as many of the defined parameters as you want, in any order, separated by commas. For example, if you want to run at `c=20` and with 300 trials in parallel instead of the default 100, you could use either of

     decisionMake(c=20, numdots=300);

or

     decisionMake(numdots=300, c=20);

both have exactly the same effect.

### List of allowed parameters and their default values:

     mu0=2, c=0, wR=0, wI=6.8, noise_sigma=0.05, deltat=0.05,
     bound=0.97, maxiter=100, numdots=80, displaytails=true

### Description of each allowed parameter:

- `mu0` constant excitation to both units

- `c`   random dots coherence

- `wR`  strength of self-excitation connection

- `wI`  strength of inhibition connection between the two units

- `noise_sigma`  standard deviation of noise added at each timestep

- `deltat`   timestep

- `bound`     value of r1 and r2 at which a decision is deemed to be emitted

- `maxiter`   Maximum number of timesteps in each simulation

- `numdots`   Number of trials to simulate in parallel

- `displaytails`  Must be one of `true` or `false`. If true, show a green line tracing each trajectory, if false show only endpoints.


