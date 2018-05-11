### Writeup - JWG2 - 5.11.18

#### Design
In constrast to assignment two I chose to work entirely outside of Jupyter notebook,
which was definitely the right call. There are many separate pieces to the Solver, and
trying to cram them all into a notebook would make it very difficult to operate.
Initially, I wanted the RHSGenerator to be integrated with the main Solver, meaning
the right hand side values would generate at runtime and have their values populated.
I wasn't able to get this to work because of the complexity of substituting
symbols with values using Sympy, so instead RHSGenerator acts independently,
and I created the dSdt equations based on its output.

Not being constrained to a notebook my reaction rate equations are in
CalcReactionRates for organization.
<br><br>
#### Difficulty
Writing the code was the easiest part! What was much harder was *understanding*
what to code.

After setting up a very simple solver it was difficult to add components such
as temperature and the density dependent reactions. Specifically, trying to make
sense of the formulas we discussed in class was challenging as I have very
little chemistry experience. Additionally, Grackle goes much more in depth with
their calculations than we do, so a large amount of the time I spent on the
project was simply reading PDFs and the Grackle source to try and figure out
what exactly I needed and what I did not.  
<br>
#### Speed and Equilibrium
As seen in [HighTemp.png](https://github.com/JakobGrootens/ChemistrySolver/blob/master/HighTemp.png) and [LowTemp.png](), having a higher temperature
creates a much more dynamic simulation with high rates of change. Conversely,
a low temperature results in relatively little change over the course of the
simulation. In both cases though, the species start in somewhat of an equilibrium,
change a great deal in the middle of the simulation, and then even out towards
an equilibrium at the end.
<br><br>
#### Interesting Regions
One region of particular interest is in the [high temperature simulation](https://github.com/JakobGrootens/ChemistrySolver/blob/master/HighTemp.png)
between T = 10^4 and 10^6. It is a moment where the simulation shifts very quickly
from a slow and steady reaction to one of enormous change, but then suddenly
almost all change stops!
