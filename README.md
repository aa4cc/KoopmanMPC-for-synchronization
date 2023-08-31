# KoopmanMPC-for-synchronization

This project demonstrates the application of Koopman MPC framework for controlled synchronization of the Frenkel-Kontorova Model published at 
[http://dx.doi.org/10.1016/j.conengprac.2023.105629](http://dx.doi.org/10.1016/j.conengprac.2023.105629)

To run the examples:
- Open the ./simulations folder in the Matlab.
- Add all folders and subfolders into the **path**.
- In the directory ./MPC_functions/, run **install_osqp.m**. This should install the OSQP solver for QP for your system. Refer to https://osqp.org/ for further information.
- The examples are in the folder ./simulation_demos/
- The comparison of MPC and KMPC are in the folder ./simulation_demos_MPC_KMPC_comparison/