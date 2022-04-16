[![Build status](https://github.com/ebranlard/wtDigiTwin-data/workflows/Tests/badge.svg)](https://github.com/ebranlard/wtDigiTwin-data/actions?query=workflow%3A%22Tests%22)


# Main Notes/Reminders

- NOTE: BldNodes in ED affects Fzi in OpenFAST and has a huge impact on Heave motion for TS

- NOTE: Cd=0. Lin models cannot be good with Cd/=0 for now

- Loads from "0" needs to be mapped at "O" using sin/cos




# Sturday 16th of April

- Issue in Loads at Origin, was that the transfer of loads needs "sin" and "cos" of theta
  TODO understand why

- Better calculation of masses using "OpenFAST" algo


- HydroSim_PrescribeLoads:
   - All OK, at 0/or O. Linear Yaw weird. Linear vel and acc show oscillations and are sensitive
     to uop, qop, and maybe amplitude of oscillations
   - Results are identical whether using "0" or "O"


# Friday 15th of April

- NOTE: Important GOTCHA: Cd, RefPointMotion, RefPointMapping, number of BldNodes, qop, uop,
    - NOTE: BldNodes in ED affects Fzi in OpenFAST and has a huge impact on Heave motion for TS
    - NOTE: Cd=0. Lin models cannot be good with Cd/=0 for now

- There was a distinciton to be made between RefPointMotion and RefPointMapping.
  hydro0 and hydroO models differ in where the mapping is done
  This affect the linearized matrix 
       Motion Mesh RefPoint was always (0,0,0)
       Mapping was always (0,0,0) suitable for "hydro0" models but for not hydroO


- HydroSim_PrescribeMotion
    - TS, Spar and SparNoRNA give perfect match for NL, Decent match for Lin as long as Cd=0!
    - TODO: Fix hydrodyn\_driver to support PRP and EDRef Motions and Loads
    - TODO (longterm): Fix lin for Cd/=0


- HydroSim_PrescribeLoads:
   - SparNoRNA: OK/Perfect for hydro0 
   - Spar: OK for hydro0
   - (DEBUG): SparNoRNA with hydroO works for 111111 but with wrong sign in Refz! That's pretty weird!
              Use 0000010 , it doesn't work well for that one

   - (DEBUG): TS Pitch does not work!!!! use pendulum_hydro.. to debug. 
                 Most likely issue with uop and qop
                 OR issue with the fact that "0" point is tilted
    > BOTH FIXED using "sin"/"cos" transfer

- HydroSim Coupled: 
    - Works well for SparNoRNA 
    - Not so well for Spar. Maybe qop? uop?
    - TODO: fix TODO above before moving to more advanced TS case
