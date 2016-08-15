# snowGlobeFluid
(Java/Processing) Fluid Globe with rb-based snowmen and particle showflakes

This repo is an eclipse project, with all libraries it needs to run (on a windows-based system). 
It should also work on other systems although it hasn't been tested.

This is a combination of a few older projects - my first graphics project where I made the original snowman animation, a java port of 
my c++ baraff-witkin-based constrained particle/rigid body sim, and a stam-inspired eulerian fluid sim with vorticity confinement.  The 
ultimate goal is to have the snowman animation continue within the snow globe until the fluid forces upon them breach some threshold
whereupon the snowmen will fly apart and be transported around along with the snowflakes as rigid bodies.  Then, once the fluid forces
drecrease below some thresholds, forces are applied to the snowmen for them to flow back together (minimize the distance from each piece
to the COM of the constituent pieces) in order to reconstitute and resume their snowball fight. 

This is a work in progress.