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
	//Sim class bindings    
	//---------------------------------------
	py::class_<Sim>(m,"Sim")
		.def(py::init< >())
		.def_readwrite("nparticles", &Sim::nparticles)
		.def_readwrite("E0", &Sim::E0)
		.def_readwrite("rho", &Sim::rho)
		.def_readwrite("nparticles", &Sim::nparticles)
		.def_readwrite("box_relative", &Sim::box_relative)
		.def_readwrite("beta", &Sim::beta)
		.def_readwrite("pressure", &Sim::pressure)
		.def_readwrite("dt", &Sim::dt)
		.def_readwrite("time", &Sim::time)
		.def_readwrite("tot_energy", &Sim::tot_energy)
		.def_readwrite("ke", &Sim::ke)
		.def_readwrite("pe", &Sim::pe)
		.def_readwrite("p", &Sim::p)
		.def_readwrite("Tinst", &Sim::Tinst)
		.def_readwrite("Pinst", &Sim::Pinst)
		.def_readwrite("Pinst_1", &Sim::Pinst_1)
		.def_readwrite("virial", &Sim::virial)
		.def_readwrite("hypervirial", &Sim::hypervirial)
		.def_readwrite("press_kin", &Sim::press_kin)
		.def_readwrite("stress", &Sim::stress)
		.def_readwrite("c_v", &Sim::c_v)
		.def_readwrite("c_p", &Sim::c_p)
		.def_readwrite("gamma_v", &Sim::gamma_v)
		.def_readwrite("beta_S", &Sim::beta_S)
		.def_readwrite("beta_T", &Sim::beta_T)
		.def_readwrite("alpha_p", &Sim::alpha_p)
		.def_readwrite("box", &Sim::box)
		.def_readwrite("box_2", &Sim::box_2)
		.def_readwrite("volume", &Sim::volume)
		.def_readwrite("integrator", &Sim::integrator)
		//add rdf properties in a bit
		//file pointer?
		.def_property("particles", &Sim::gparticles, &Sim::sparticles)
		.def_property("potential", &Sim::gpotential, &Sim::spotential)
		.def_property("thermostat", &Sim::gthermostat, &Sim::sthermostat)
		.def_property("barostat", &Sim::gbarostat, &Sim::sbarostat)
		//methods
		.def("forces", &Sim::forces)
		.def("potential_energy", &Sim::potential_energy)
		.def("kinetic_energy", &Sim::kinetic_energy)
		.def("total_energy", &Sim::total_energy)
		.def("total_momentum", &Sim::total_momentum)
		.def("rescale_velocities", &Sim::rescale_velocities)
		.def("image_distance", &Sim::image_distance)
		.def("remap", &Sim::remap)
		.def("init", &Sim::init)
		.def("md_step", &Sim::md_step)
		.def("md_verlet", &Sim::md_verlet)
		.def("md_langevin", &Sim::md_langevin)
		.def("md_andersen", &Sim::md_andersen)
		.def("md_nosehooverlangevin_NVT", &Sim::md_nosehooverlangevin_NVT)
		.def("md_andersen_NPH", &Sim::md_andersen_NPH)
		.def("md_andersen_stochastic_NPT", &Sim::md_andersen_stochastic_NPT)
		.def("md_andersen_stochastic_nhlthermo_NPT", &Sim::md_andersen_stochastic_nhlthermo_NPT)
		.def("langevin_thermo", &Sim::langevin_thermo)
		.def("langevin_baro", &Sim::langevin_baro)
		.def("propagate_momenta_half", &Sim::propagate_momenta_half)
		.def("propagate_position_half", &Sim::propagate_position_half)
		.def("propagate_momenta_xi", &Sim::propagate_momenta_xi)
		.def("langevin_xi", &Sim::langevin_xi)
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
		.def_readwrite("fi", &Potential::fi)
		.def_readwrite("fj", &Potential::fj)
		.def("forces", &Potential::forces)
		.def("potential_energy", &Potential::potential_energy)
		.def("init", &Potential::init)
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
		.def_readwrite("dElangevin", &Thermostat::dElangevin)
		.def("init", &Thermostat::init)
		.def("init_nhl", &Thermostat::init_nhl)
		;

	py::class_<Barostat>(m, "Barostat")
		.def(py::init< >())
		.def_readwrite("isotropic", &Barostat::isotropic)
		.def_readwrite("pv", &Barostat::pv)
		.def_readwrite("pmass", &Barostat::pmass)
		.def_readwrite("lgamma", &Barostat::lgamma)
		.def_readwrite("lc1", &Barostat::lc1)
		.def_readwrite("lc2", &Barostat::lc2)
		.def("init", &Barostat::init)
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