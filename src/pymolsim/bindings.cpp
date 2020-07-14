/*
Python bindings for LJ MD code
*/

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <vector>
#include <string>
#include "header.h"

namespace py = pybind11;
using namespace std;

PYBIND11_MODULE(cmolsim, m) {
    py::options options;
    options.disable_function_signatures();

	//---------------------------------------
	//System class bindings    
	//---------------------------------------
	py::class_<System>(m,"System")
		.def(py::init< >())
		.def_readwrite("ncycle1", &System::ncycle1)
		.def_readwrite("ncycle2", &System::ncycle2)
		.def_readwrite("nparticles", &System::nparticles)
		.def_readwrite("E0", &System::E0)
		.def_readwrite("rho", &System::rho)
		.def_readwrite("nparticles", &System::nparticles)
		.def_readwrite("box_relative", &System::box_relative)
		.def_readwrite("beta", &System::beta)
		.def_readwrite("pressure", &System::pressure)
		.def_readwrite("dt", &System::dt)
		.def_readwrite("time", &System::time)
		.def_readwrite("tot_energy", &System::tot_energy)
		.def_readwrite("dElangevin", &System::dElangevin)
		.def_readwrite("Tinst", &System::Tinst)
		.def_readwrite("Pinst", &System::Pinst)
		.def_readwrite("Pinst_1", &System::Pinst_1)
		.def_readwrite("virial", &System::virial)
		.def_readwrite("hypervirial", &System::hypervirial)
		.def_readwrite("press_kin", &System::press_kin)
		.def_readwrite("stress", &System::stress)
		.def_readwrite("c_v", &System::c_v)
		.def_readwrite("c_p", &System::c_p)
		.def_readwrite("gamma_v", &System::gamma_v)
		.def_readwrite("beta_S", &System::beta_S)
		.def_readwrite("beta_T", &System::beta_T)
		.def_readwrite("alpha_p", &System::alpha_p)
		.def_readwrite("box", &System::box)
		.def_readwrite("box_2", &System::box_2)
		.def_readwrite("volume", &System::volume)
		.def_readwrite("integrator", &System::integrator)
		//add rdf properties in a bit
		//file pointer?
		;







#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}