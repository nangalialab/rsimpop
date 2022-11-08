# CHANGES
## 2.2.3
Clean up of vignette documentation
## 2.2.2
Now allow the use of a continuous driver acquisition process in conjunction with a time varying target population size trajectory
## 2.2.1
Optimisation of multiple driver handling: 
* Much faster run_driver_process_sim when drivers_per_year becomes high (>10)  and hopefully similar performance otherwise
* Each added driver now has a unique driver ID that is preserved between runs of run_driver_process_sim and continue_driver_process_sim
* run_driver_process_sim now can take either a function for generating selection coefficients or an explicit vector of coefficients.  If you provide a vector then make sure it is at least 200,000 elements long.
* plot_tree_events has an additional parameter that will enable consistent colour/shape plotting of drivers between two runs on the same individual

## 2.2.0

Changes contributed by Karim Mane

* Added unique identifier (uid) column to the events table. 
* combine_simpop is now redundant.

## 2.1.2

* run_neutral_trajectory now works when initialised with a non-null simpop object.

## 2.1.1

* Added missing 2.1.0  changes to CHANGES.md

## 2.1.0

* Removal of redundant tree plotting functions

## 2.0.4

* Switched from multiplicative to additive selection model for compound drivers
* Added parameter to plot_tree_events to show selection coefficient on legend.

## 2.0.3

* Fixed bug in enhanced handling of drivers.

## 2.0.2

* Enhanced handling of drivers so that there is no set cap on the number of drivers.

## 2.0.1

* Initial Release
