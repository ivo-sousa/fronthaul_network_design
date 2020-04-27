# Fronthaul Design for Wireless Networks: a Simulation Tool
This is the location of a simulator, developed using the MATLAB programming language, which has the following goals: 
(i) to test and compare the performance of Microwave Radio Transmission (MRT), Free Space Optics (FSO) and Fiber Optics (FO) technologies under user-specified equipment characteristics and link conditions, thus enabling to determine the most cost-effective solution to connect two points {Link Design Algorithm}; 
(ii) to optimally find a fronthaul topology for wireless networks, given the Remote Radio Heads (RRHs) location — namely, the required number of Baseband Units (BBUs) and where they should be positioned in order to minimize the overall network costs {Network Planning Algorithm}.

The algorithms associated with these goals, Link Design Algorithm and Network Planning Algorithm, can be found in the files 'link_design_algorithm.m' and 'network_planning_algorithm.m', respectively, within the 'Simulator' folder.

For more information, namely detailed descriptions of the aforementioned algorithms and some examples of results that can be obtained using the simulator, please refer to the 'Algorithms Details & Examples.pdf' file.

## License
This simulator is released under MIT License. See the bundled [LICENSE](LICENSE) file for details.
