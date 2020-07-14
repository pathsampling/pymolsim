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

	py::class_<Potential>(m, "Potential")
		.def(py::init< >())
		.def_readwrite("sigma", &Potential::sigma)
		.def_readwrite("epsilon", &Potential::epsilon)
		.def_readwrite("c1", &Potential::c1)
		.def_readwrite("c2", &Potential::c2)
		.def_readwrite("c3", &Potential::c3)
		.def_readwrite("c4", &Potential::c4)
		.def_readwrite("c5", &Potential::c5)
		.def_readwrite("rmin2", &Potential::rmin2)
		.def_readwrite("rmax2", &Potential::rmax2)
		.def_readwrite("virial", &Potential::virial)
		.def_readwrite("hypervirial", &Potential::hypervirial)
		.def_readwrite("stress", &Potential::stress)
		;

	py::class_<Particle>(m, "Particle")
		.def(py::init< >())
		.def_readwrite("r", &Particle::r)
		.def_readwrite("v", &Particle::v)
		.def_readwrite("f", &Particle::f)
		.def_readwrite("energy", &Particle::energy)
		.def_readwrite("mass", &Particle::mass)
		.def_readwrite("c", &Particle::c)
		.def_readwrite("type", &Particle::type)
		.def_readwrite("c_z", &Particle::c_z)
		.def_readwrite("phi", &Particle::phi)
		;

	py::class_<Thermostat>(m, "Thermostat")
		.def(py::init< >())
		.def_readwrite("lgamma", &Thermostat::lgamma)
		.def_readwrite("lc1", &Thermostat::lc1)
		.def_readwrite("lc2", &Thermostat::lc2)
		.def_readwrite("anu", &Thermostat::anu)
		.def_readwrite("nhlgamma", &Thermostat::nhlgamma)
		.def_readwrite("nhlmu", &Thermostat::nhlmu)
		.def_readwrite("nhlc1", &Thermostat::nhlc1)
		.def_readwrite("nhlc2", &Thermostat::nhlc2)
		.def_readwrite("nhlxi", &Thermostat::nhlxi)
		;

	py::class_<Barostat>(m, "Barostat")
		.def(py::init< >())
		.def_readwrite("isotropic", &Barostat::isotropic)
		.def_readwrite("pv", &Barostat::pv)
		.def_readwrite("pmass", &Barostat::pmass)
		.def_readwrite("lgamma", &Barostat::lgamma)
		.def_readwrite("lc1", &Barostat::lc1)
		.def_readwrite("lc2", &Barostat::lc2)
		;

	//read input will be deprecated
	//read init_structure and read_initlammpsdump
    //needs binding






#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}