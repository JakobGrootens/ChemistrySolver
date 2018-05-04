### Progress Report - JWG2 - 5.4.18

So far I've got most of the structure down, including a Right-Hand-Side generator, and many of the reaction rates from Grackle. The results of a preliminary run are in [simulation.png](https://github.com/JakobGrootens/ChemistrySolver/blob/master/simulation.png)

I have the 19 main reactions from Grackle source, but still need to implement four more. Two of those I simply haven't found yet (HI + HI) and (HI + HeI), and the other two are density dependent. For the density dependent reactions I also need to use the ideal gas law formula from lecture:

e = \frac{kT}{(\gamma - 1) \mu m\_{\rm H}} 

My only concerns at this point are being able to test the "correctness" of my code, and updating my electron count properly. I'm not sure how to calculate de for the integrator since the recommended approach is to recalculate n_e each iteration. My current intuition is to manually change e in the integrator's vector so I can graph it properly, or remove e from the state vector entirely.
