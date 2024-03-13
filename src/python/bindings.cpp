#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <iostream>
#include <ctime>
#include <random>
#include <vector>
#include <chrono>
#include <cmath>
#include <any>
#include "dd/MDDPackage.hpp"



namespace py = pybind11;


using Instruction = std::tuple< std::string, bool, std::vector<int>, std::string, std::vector<int>, py::object, std::tuple<std::vector<dd::QuantumRegister>, std::vector<dd::Control::Type>> >;
using Circuit = std::vector<Instruction>;
using Circuit_info = std::tuple<unsigned int, std::vector<size_t>, Circuit>;

using Noise = std::tuple<double, double>;
using NoiseType = std::map< std::variant<std::string, std::vector<int>>, Noise>;
using NoiseModel = std::map<std::string, NoiseType>;
using CVec = std::vector<std::complex<double>>;

// =======================================================================================================
// =======================================================================================================
// =======================================================================================================
// ============@@@@@@@@@@@@@@@@@@@@ PRINTING FUNCTION @@@@@@@@@@@@@@@@@@@@================================
// =======================================================================================================
// =======================================================================================================


// Function to print NoiseType
void printNoiseType(const NoiseType& noiseType, int indent = 0) {
    for (const auto& [key, noise] : noiseType) {
        // Print the key
        for (int i = 0; i < indent; ++i)
            std::cout << "    ";
        if (std::holds_alternative<std::string>(key))
            std::cout << std::get<std::string>(key) << " -> ";
        else {
            std::cout << "[";
            for (int val : std::get<std::vector<int>>(key))
                std::cout << val << " ";
            std::cout << "] -> ";
        }

        // Print the noise tuple
        std::cout << "(" << std::get<0>(noise) << ", " << std::get<1>(noise) << ")" << std::endl;
    }
}

// Function to print NoiseModel
void printNoiseModel(const NoiseModel& noiseModel) {
    for (const auto& [key, noiseType] : noiseModel) {
        std::cout << key << ":" << std::endl;
        printNoiseType(noiseType, 1);
    }
}

void printCvec(CVec vector){
    for (const auto& cn: vector) {
        std::cout << cn << std::endl;
    }
}


void printCircuit(const Circuit& circuit) {
    for (const auto& instruction : circuit) {
        auto [tag, dag, dims, gate_type, target_qudits, params, control_set] = instruction;
        std::cout << "Tag: " << tag << std::endl;
        std::cout << "Dag: " << dag << std::endl;
        std::cout << "Dimensions: ";
        for (const auto& dim : dims) {
            std::cout << dim << " ";
        }
        std::cout << std::endl;
        std::cout << "Gate Type: " << gate_type << std::endl;
        std::cout << "Target Qudits: ";
        for (const auto& qubit : target_qudits) {
            std::cout << qubit << " ";
        }
        std::cout << std::endl;

        // Printing control_set
        auto [control1, control2] = control_set;
        std::cout << "Control Set: ";
        for (const auto& control : control1) {
            std::cout << control << " ";
        }
        std::cout << "| ";
        for (const auto& control : control2) {
            std::cout << control << " ";
        }
        std::cout << std::endl;
    }
}

// =======================================================================================================
// =======================================================================================================
// =======================================================================================================
// ============================================= HELPER FUNCTION =========================================
// =======================================================================================================
// =======================================================================================================

bool is_none_or_empty(const py::object &obj) {
    if (obj.is_none())
        return true;

    if (py::isinstance<py::sequence>(obj)) {
        auto seq = obj.cast<py::sequence>();
        return seq.size() == 0;
    }

    // Add additional checks for other types if needed

    return false;
}

bool checkDim(const std::vector<int>& dims, const std::variant<int, std::vector<int>>& target) {
    if (std::holds_alternative<int>(target)) {
        // If target is a single integer
        if (dims.size() != 1) {
            return false; // Different sizes, not exactly equal
        }
        return dims[0] == std::get<int>(target);
    } else {
        // If target is a vector
        const auto& targetVec = std::get<std::vector<int>>(target);
        if (dims.size() != targetVec.size()) {
            return false; // Different sizes, not exactly equal
        }
        for (size_t i = 0; i < dims.size(); ++i) {
            if (dims[i] != targetVec[i]) {
                return false; // Different elements, not exactly equal
            }
        }
        return true; // All elements are the same, exactly equal
    }
}

// Function to convert C++ vector of complex numbers to a Python list
py::list complex_vector_to_list(const CVec& vec) {
    py::list pyList;
    for (const auto& elem : vec) {
        try {
            // Convert std::complex<double> to Python complex object
            py::object pyComplex = py::cast(elem);
            // Append Python complex object to the list
            pyList.append(pyComplex);
        } catch (const std::exception& e) {
            // Handle any exceptions
            std::cerr << "Error appending Python object: " << e.what() << std::endl;
        }
    }
    return pyList;
}

// =======================================================================================================
// =======================================================================================================
// =======================================================================================================
// ===============::::::::::::::::::::::PARSIN FUNCTIONS ::::::::::::::::::::::===========================
// =======================================================================================================
// =======================================================================================================



Circuit_info readCircuit(py::object &circ) {
    Circuit result;
    std::cout << "reading" << std::endl;

    unsigned int num_qudits = circ.attr("_num_qudits").cast<unsigned int>();
    std::vector<size_t> dimensions = circ.attr("_dimensions").cast<std::vector<size_t>>();
    
    std::cout << "No need to store this string";
    
    bool result_empty = py::isinstance<py::none>(circ.attr("instructions"));
    std::cout << "Is empty: " << std::boolalpha << result_empty << std::endl;
    
    // Get Python iterable
    py::iterator it = py::iter(circ.attr("instructions"));
    
    // Iterate over the Python iterable
    while (it != py::iterator::sentinel()) {
        // Get the current Python object
        //py::object obj = *it;
        py::handle obj_handle = *it;
        py::object obj = py::reinterpret_borrow<py::object>(obj_handle);
       
        // Print the string representation
        std::string tag = obj.attr("qasm_tag").cast<std::string>();
        std::cout << tag << std::endl;

        bool dagger = obj.attr("dagger").cast<bool>();

	  // Get the string representation of the Python enum value
	  std::string gate_type = py::cast<std::string>(obj.attr("gate_type").attr("name"));

        //std::string gate_type = obj.attr("gate_type").cast<std::string>();
        std::cout << gate_type << std::endl;

        // Extracting dimensions
        py::object dims_obj = obj.attr("_dimensions");
        std::vector<int> dims;
        if (py::isinstance<py::int_>(dims_obj)) {
            dims.push_back(py::cast<int>(dims_obj));
        } else if (py::isinstance<py::list>(dims_obj)) {
            dims = py::cast<std::vector<int>>(dims_obj);
        }

        // Extracting target_qudits
        py::object target_qudits_obj = obj.attr("_target_qudits");
        std::vector<int> target_qudits;
        if (py::isinstance<py::int_>(target_qudits_obj)) {
            target_qudits.push_back(py::cast<int>(target_qudits_obj));
        } else if (py::isinstance<py::list>(target_qudits_obj)) {
            target_qudits = py::cast<std::vector<int>>(target_qudits_obj);
        }
        
        std::cout << "vector ints"<< std::endl;
        
        py::object params;
        if (is_none_or_empty(obj.attr("_params"))){
        	std::cout << "params true"<< std::endl;
        }
        else{
        std::cout << "params false"<< std::endl;
        	params = obj.attr("_params");
        }
        std::cout << "params"<< std::endl;
        
        std::tuple<std::vector<dd::QuantumRegister>, std::vector<dd::Control::Type>> control_set = {};
        if (is_none_or_empty(obj.attr("_controls_data"))){
        	        
        	        std::cout << "control true"<< std::endl;
        }
        else{
            std::cout << "control false"<< std::endl;
            py::object controls_data = obj.attr("_controls_data");
	        std::vector<dd::QuantumRegister> indices = controls_data.attr("indices").cast<std::vector<dd::QuantumRegister>>();
	        std::vector<dd::Control::Type> ctrl_states = controls_data.attr("ctrl_states").cast<std::vector<dd::Control::Type>>();
	         control_set = std::make_tuple(indices, ctrl_states);
        }
        
        result.push_back(std::make_tuple(tag, dagger, dims, gate_type, target_qudits, params, control_set));
        
    // Increment the iterator
        ++it;
    }
    
	return std::make_tuple(num_qudits, dimensions, result);
    
}

NoiseModel parse_noise_model(const py::dict& noise_model) {
    NoiseModel newNoiseModel;
    std::cout << "Noise Model"<< std::endl;
    
    for (const auto& gate : noise_model) {
        auto gateName = gate.first.cast<std::string>();
        std::cout << gateName<< std::endl;
        
        py::dict gateNoise = gate.second.cast<py::dict>();
        std::cout << gateNoise << std::endl;
        
        NoiseType newNoiseType;
        std::cout << "Noise Type"<< std::endl;
        
        std::variant<std::string, std::vector<int>> noiseSpread; // Declared outside the if blocks
        
        for (const auto& noiseTypesPair : gateNoise) {
            if (py::isinstance<py::str>(noiseTypesPair.first)) {
                std::cout << "inside"<< std::endl;
                std::string noiseSpreadString = noiseTypesPair.first.cast<std::string>();
                noiseSpread = noiseSpreadString;
                // Handle string case
            } else if (py::isinstance<py::tuple>(noiseTypesPair.first)) {
                std::cout << "outside"<< std::endl;
                std::vector<int> noiseSpreadTuple = noiseTypesPair.first.cast<std::vector<int>>();
                noiseSpread = noiseSpreadTuple;
            }
                    std::cout << "noiseeeeee"<< std::endl;
                    
                double depo = noiseTypesPair.second.attr("probability_depolarizing").cast<double>();
	       double deph = noiseTypesPair.second.attr("probability_dephasing").cast<double>();
	       std::tuple<double, double> noiseProb = std::make_tuple(depo, deph);
            //std::tuple<double, double> noiseProb = noiseTypesPair.second.cast<std::tuple<double,double>>();
            newNoiseType[noiseSpread] = noiseProb;
            std::cout << "Noise Spread"<< std::endl;
        }
        newNoiseModel[gateName] = newNoiseType;
    }
    return newNoiseModel;
}

Circuit generateCircuit(const Circuit_info& circuitInfo, const NoiseModel& noiseModel) {
    // Get current time in milliseconds
    auto currentTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    const auto& [num_qudits, dimensions, circuit] = circuitInfo;
    // Calculate seed
    std::random_device rd;
    std::mt19937_64 gen(rd() + currentTimeMs);

    Circuit noisyCircuit;

    for (const Instruction& instruction : circuit) {
        noisyCircuit.push_back(instruction);

        const auto& [tag, dag, dims, gate_type, target_qudits, params, control_set] = instruction;
        std::vector<int> referenceLines(target_qudits.begin(), target_qudits.end());

        // If controls is not empty, add its first field to reference_lines
        if (!(std::get<0>(control_set).empty()) && (std::get<1>(control_set).empty()) ) {
            auto [ctrl_dits, levels] = control_set; // Decompose the tuple

            // Choose which vector to append to reference_lines based on your logic
            // For example, appending control_vec1 here
            referenceLines.insert(referenceLines.end(), ctrl_dits.begin(), ctrl_dits.end());
        }

        if (noiseModel.find(tag) != noiseModel.end()) {
            for (const auto& mode_noise : noiseModel.at(tag)) {
                auto mode = mode_noise.first;
                auto noise_info = mode_noise.second;
                double depo = std::get<0>(noise_info);
                double deph = std::get<1>(noise_info);

                std::discrete_distribution<int> x_dist({1.0 - depo, depo});
                std::discrete_distribution<int> z_dist({1.0 - deph, deph});

                int x_choice = x_dist(gen);
                int z_choice = z_dist(gen);
                std::cout << "X CHOICE  "<< x_choice<< std::endl;
                std::cout << "Z CHOICE  "<< z_choice<< std::endl;

                if (x_choice == 1 || z_choice == 1) {
                    std::cout << "DAYUUM  "<< std::endl;

                    std::vector<int> qudits;

                    if (std::holds_alternative<std::vector<int>>(mode)) {
                        std::cout << "qudits is int vec  "<< std::endl;
                        qudits = std::get<std::vector<int>>(mode);

                    } else if (std::holds_alternative<std::string>(mode)) {
                        std::cout << "string mode "<< std::endl;
                        std::string mode_str = std::get<std::string>(mode);
                        if (mode_str == "local") {
                            qudits = referenceLines;
                        } else if (mode_str == "all") {
                            for (int i = 0; i < num_qudits; ++i)
                                qudits.push_back(i);
                        } else if (mode_str == "nonlocal") {
                            assert(gate_type == "TWO" || gate_type == "MULTI");
                            qudits = referenceLines;
                        } else if (mode_str == "control") {
                            assert(gate_type == "TWO");
                            qudits.push_back(target_qudits.at(0));
                        } else if (mode_str == "target") {
                            assert(gate_type == "TWO");
                            qudits.push_back(target_qudits.at(1));
                        }
                    }

                    if (x_choice == 1) {
                        for (auto dit : qudits) {
                            if (tag == "rxy" || tag == "rz" || tag == "virtrz") {
                                py::list params_new;

                                // Retrieve field 0 and 1 from params
                                py::object field_0 = params.attr("__getitem__")(0);
                                py::object field_1 = params.attr("__getitem__")(1);

                                // Convert field_0 and field_1 to integers (assuming they are integers)
                                int value_0 = py::cast<int>(field_0);
                                int value_1 = py::cast<int>(field_1);

                                // Create a new list and append value_0 and value_1
                                params_new.append(value_0);
                                params_new.append(value_1);

                                // Append pi and pi/2
                                double pi = 3.14159265358979323846;
                                double pi_over_2 = pi / 2.0;
                                params_new.append(py::float_(pi));
                                params_new.append(py::float_(pi_over_2));
                                Instruction new_inst = std::make_tuple("rxy", false, dims, "SINGLE", std::vector<int>{dit}, py::cast<py::object>(params_new), std::tuple<std::vector<dd::QuantumRegister>, std::vector<dd::Control::Type>>());
                                noisyCircuit.push_back(new_inst);
                            } else {
                                py::object params_new;

                                // Construct params_new for other gates (e.g., x)
                                // Example: params_new = some_value;
                                Instruction new_inst = std::make_tuple("x", false, dims, "SINGLE", std::vector<int>{dit}, params_new, std::tuple<std::vector<dd::QuantumRegister>, std::vector<dd::Control::Type>>());
                                noisyCircuit.push_back(new_inst);
                            }
                        }
                    }

                    if (z_choice == 1) {
                        for (auto dit : qudits) {
                            if (tag == "rxy" || tag == "rz" || tag == "virtrz") {
                                py::list params_new;

                                // Retrieve field 0 and 1 from params
                                py::object field_0 = params.attr("__getitem__")(0);
                                py::object field_1 = params.attr("__getitem__")(1);

                                // Convert field_0 and field_1 to integers (assuming they are integers)
                                int value_0 = py::cast<int>(field_0);
                                int value_1 = py::cast<int>(field_1);

                                // Create a new list and append value_0 and value_1
                                params_new.append(value_0);
                                params_new.append(value_1);

                                // Append pi and pi/2
                                double pi = 3.14159265358979323846;
                                params_new.append(py::float_(pi));
                                Instruction newInst = std::make_tuple("rz", false, dims, "SINGLE", std::vector<int>{dit}, py::cast<py::object>(params_new), std::tuple<std::vector<dd::QuantumRegister>, std::vector<dd::Control::Type>>());
                                noisyCircuit.push_back(newInst);
                            } else {
                                py::object paramsNew;

                                // Construct params_new for other gates (e.g., z)
                                // Example: params_new = some_value;
                                Instruction newInst = std::make_tuple("z", false, dims, "SINGLE", std::vector<int>{dit}, paramsNew, std::tuple<std::vector<dd::QuantumRegister>, std::vector<dd::Control::Type>>());
                                noisyCircuit.push_back(newInst);
                            }
                        }
                    }
                }
            }
        }
    }

    return noisyCircuit;
}



// =======================================================================================================
// =======================================================================================================
// =======================================================================================================
// ===============()()()()()()()()()()() SIMULATION FUNCTIONS ()()()()()()()()()==========================
// =======================================================================================================
// =======================================================================================================
// =======================================================================================================
// =======================================================================================================
/*
 * SUPPORTED GATES AT THE MOMENT UNTIL DIMENSION 5
"csum": "csum",
"cx": "cx",
"h": "h",
"rxy": "r",
"rz": "rz",
"virtrz": "virtrz",
"s": "s",
"x": "x",
"z": "z"
 */
using ddpkg = std::unique_ptr<dd::MDDPackage>;

dd::MDDPackage::mEdge getGate(const ddpkg& dd, const Instruction& instruction){
    std::cout << "getting your op "<< std::endl;
	const auto& [tag, dag, dims, gate_type, target_qudits, params, control_set] = instruction;
    dd::MDDPackage::mEdge gate;
    auto numberRegs= static_cast<dd::QuantumRegisterCount>(dd->numberOfQuantumRegisters);

    dd::QuantumRegister tq = 0;
    if(gate_type=="SINGLE"){
         tq = static_cast<dd::QuantumRegister>(target_qudits.at(0));
    }

    dd::Controls controlSet{};
    if (!(std::get<0>(control_set).empty()) && (std::get<1>(control_set).empty()) ) {
        std::vector<dd::QuantumRegister> ctrlQudits = std::get<0>(control_set);
        std::vector<dd::Control::Type> ctrlLevels = std::get<1>(control_set);
        for( uint i=0; i < ctrlQudits.size(); i++){
            const dd::Control c{ctrlQudits.at(i), ctrlLevels.at(i)};
            controlSet.insert(c);
        }
    }
    std::cout << "selecting   "<< std::endl;
    if (tag == "rxy") {
        // Handle rxy tag with dimension 2
        auto pl = params.cast<py::list>();
        auto leva = pl[0].cast<size_t>();
        auto levb = pl[1].cast<size_t>();
        auto theta = pl[2].cast<double>();
        auto phi = pl[3].cast<double>();

        if (checkDim(dims, 2)) {
            dd::GateMatrix matrix = dd::RXY(theta, phi);
            gate = dd->makeGateDD<dd::GateMatrix>(matrix, numberRegs, controlSet, tq);

        } else if (checkDim(dims, 3)) {
            // Handle rxy tag with dimension 3
            dd::TritMatrix matrix = dd::RXY3(theta, phi, leva, levb);
            gate = dd->makeGateDD<dd::TritMatrix>(matrix, numberRegs, controlSet, tq);

        } else if (checkDim(dims, 4)) {
            // Handle rxy tag with dimension 4
            dd::QuartMatrix matrix = dd::RXY4(theta, phi, leva, levb);
            gate = dd->makeGateDD<dd::QuartMatrix>(matrix, numberRegs, controlSet, tq);

        } else if (checkDim(dims, 5)) {
            // Handle rxy tag with dimension 5
            dd::QuintMatrix matrix = dd::RXY5(theta, phi, leva, levb);
            gate = dd->makeGateDD<dd::QuintMatrix>(matrix, numberRegs, controlSet, tq);
        }
    } else if (tag == "rz") {
        auto pl = params.cast<py::list>();

        auto leva = pl[0].cast<size_t>();
        auto levb = pl[1].cast<size_t>();
        auto phi = pl[2].cast<double>();

        if (checkDim(dims, 2)) {
            dd::GateMatrix matrix = dd::RZ(phi);
            gate = dd->makeGateDD<dd::GateMatrix>(matrix, numberRegs, controlSet, tq);

        } else if (checkDim(dims, 3)) {
            // Handle rxy tag with dimension 3
            dd::TritMatrix matrix = dd::RZ3(phi, leva, levb);
            gate = dd->makeGateDD<dd::TritMatrix>(matrix, numberRegs, controlSet, tq);

        } else if (checkDim(dims, 4)) {
            // Handle rxy tag with dimension 4
            dd::QuartMatrix matrix = dd::RZ4(phi, leva, levb);
            gate = dd->makeGateDD<dd::QuartMatrix>(matrix, numberRegs, controlSet, tq);

        } else if (checkDim(dims, 5)) {
            // Handle rxy tag with dimension 5
            dd::QuintMatrix matrix = dd::RZ5(phi, leva, levb);
            gate = dd->makeGateDD<dd::QuintMatrix>(matrix, numberRegs, controlSet, tq);
        }
    } else if (tag == "virtrz") {
        auto pl = params.cast<py::list>();

        auto leva = pl[0].cast<size_t>();
        auto phi = pl[1].cast<double>();

        if (checkDim(dims, 2)) {
            dd::GateMatrix matrix = dd::VirtRZ(phi, leva);
            gate = dd->makeGateDD<dd::GateMatrix>(matrix, numberRegs, controlSet, tq);

        } else if (checkDim(dims, 3)) {
            // Handle rxy tag with dimension 3
            dd::TritMatrix matrix = dd::VirtRZ3(phi, leva);
            gate = dd->makeGateDD<dd::TritMatrix>(matrix, numberRegs, controlSet, tq);

        } else if (checkDim(dims, 4)) {
            // Handle rxy tag with dimension 4
            dd::QuartMatrix matrix = dd::VirtRZ4(phi, leva);
            gate = dd->makeGateDD<dd::QuartMatrix>(matrix, numberRegs, controlSet, tq);

        } else if (checkDim(dims, 5)) {
            // Handle rxy tag with dimension 5
            dd::QuintMatrix matrix = dd::VirtRZ5(phi, leva);
            gate = dd->makeGateDD<dd::QuintMatrix>(matrix, numberRegs, controlSet, tq);
        }
    } else if (tag == "x") {
        if (checkDim(dims, 2)) {
            dd::GateMatrix matrix = dd::Xmat;
            gate = dd->makeGateDD<dd::GateMatrix>(matrix, numberRegs, controlSet, tq);

        } else if (checkDim(dims, 3)) {
            // Handle rxy tag with dimension 3
            dd::TritMatrix matrix = dd::X3;
            gate = dd->makeGateDD<dd::TritMatrix>(matrix, numberRegs, controlSet, tq);

        } else if (checkDim(dims, 4)) {
            // Handle rxy tag with dimension 4
            dd::QuartMatrix matrix = dd::X4;
            gate = dd->makeGateDD<dd::QuartMatrix>(matrix, numberRegs, controlSet, tq);

        } else if (checkDim(dims, 5)) {
            // Handle rxy tag with dimension 5
            dd::QuintMatrix matrix = dd::X5;
            gate = dd->makeGateDD<dd::QuintMatrix>(matrix, numberRegs, controlSet, tq);
        }
    } else if (tag == "s") {
        if (checkDim(dims, 2)) {
            dd::GateMatrix matrix = dd::Smat;
            gate = dd->makeGateDD<dd::GateMatrix>(matrix, numberRegs, controlSet, tq);

        } else if (checkDim(dims, 3)) {
            // Handle rxy tag with dimension 3
            dd::TritMatrix matrix = dd::S3();
            gate = dd->makeGateDD<dd::TritMatrix>(matrix, numberRegs, controlSet, tq);

        } else if (checkDim(dims, 4)) {
            // Handle rxy tag with dimension 4
            dd::QuartMatrix matrix = dd::S4();
            gate = dd->makeGateDD<dd::QuartMatrix>(matrix, numberRegs, controlSet, tq);

        } else if (checkDim(dims, 5)) {
            // Handle rxy tag with dimension 5
            dd::QuintMatrix matrix = dd::S5();
            gate = dd->makeGateDD<dd::QuintMatrix>(matrix, numberRegs, controlSet, tq);
        }
    } else if (tag == "z") {
        if (checkDim(dims, 2)) {
            dd::GateMatrix matrix = dd::Zmat;
            gate = dd->makeGateDD<dd::GateMatrix>(matrix, numberRegs, controlSet, tq);

        } else if (checkDim(dims, 3)) {
            // Handle rxy tag with dimension 3
            dd::TritMatrix matrix = dd::Z3();
            gate = dd->makeGateDD<dd::TritMatrix>(matrix, numberRegs, controlSet, tq);

        } else if (checkDim(dims, 4)) {
            // Handle rxy tag with dimension 4
            dd::QuartMatrix matrix = dd::Z4();
            gate = dd->makeGateDD<dd::QuartMatrix>(matrix, numberRegs, controlSet, tq);

        } else if (checkDim(dims, 5)) {
            // Handle rxy tag with dimension 5
            dd::QuintMatrix matrix = dd::Z5();
            gate = dd->makeGateDD<dd::QuintMatrix>(matrix, numberRegs, controlSet, tq);
        }
    } else if (tag == "h") {
        if (checkDim(dims, 2)) {
            dd::GateMatrix matrix = dd::H();
            gate = dd->makeGateDD<dd::GateMatrix>(matrix, numberRegs, controlSet, tq);

        } else if (checkDim(dims, 3)) {
            // Handle rxy tag with dimension 3
            dd::TritMatrix matrix = dd::H3();
            gate = dd->makeGateDD<dd::TritMatrix>(matrix, numberRegs, controlSet, tq);

        } else if (checkDim(dims, 4)) {
            // Handle rxy tag with dimension 4
            dd::QuartMatrix matrix = dd::H4();
            gate = dd->makeGateDD<dd::QuartMatrix>(matrix, numberRegs, controlSet, tq);

        } else if (checkDim(dims, 5)) {
            // Handle rxy tag with dimension 5
            dd::QuintMatrix matrix = dd::H5();
            gate = dd->makeGateDD<dd::QuintMatrix>(matrix, numberRegs, controlSet, tq);
        }
    } else if (tag == "cx") {
        std::cout << "IN cex "<< std::endl;
        auto pl = params.cast<py::list>();

        auto leva = pl[0].cast<size_t>();
        auto levb = pl[1].cast<size_t>();
        auto ctrlLev = pl[2].cast<dd::Control::Type>();
        auto phi = pl[3].cast<dd::fp>();

        auto cReg = static_cast<dd::QuantumRegister>(target_qudits.at(0));
        auto target = static_cast<dd::QuantumRegister>(target_qudits.at(1));
        std::cout << "about to cex "<< std::endl;
        return dd->CEX(numberRegs, ctrlLev, phi, leva, levb, cReg, target, dag);
    } else if (tag == "csum") {
        //CSUM(QuantumRegisterCount n, QuantumRegister cReg, QuantumRegister target, bool isDagger = false)
        auto cReg = static_cast<dd::QuantumRegister>(target_qudits.at(0));
        auto target = static_cast<dd::QuantumRegister>(target_qudits.at(1));
        return dd->CSUM(numberRegs, cReg, target, dag);
    }
    if(dag){
        gate = dd->conjugateTranspose(gate);
    }
    return gate;
}



CVec ddsimulator(dd::QuantumRegisterCount numLines, const std::vector<size_t>& dims, const Circuit& circuit){
    std::cout << "ddsim 0  "<< std::endl;
    std::cout << "NUM LINES "<< numLines << std::endl;
    std::cout << "dims " << std::endl;
    for (size_t i = 0; i < dims.size(); ++i) {
        std::cout << dims[i];
        if (i != dims.size() - 1) {
            std::cout << ", ";
        }
    }

    const ddpkg dd = std::make_unique<dd::MDDPackage>(numLines, dims);
    auto psi = dd->makeZeroState(numLines);
    std::cout << "going to loop "<< std::endl;
    for (const Instruction& instruction : circuit) {
        dd::MDDPackage::mEdge gate = getGate(dd, instruction);
        dd->getVectorizedMatrix(gate);
        std::cout << "EPPAAAAAAAAAAAA "<< std::endl;
	    psi = dd->multiply(gate, psi);
          
    }
    return dd->getVector(psi);
}


py::list stateVectorSimulation(py::object &circ, py::object & noiseModel){
//py::list stateVectorSimulation(){
    std::cout << "INSIDE" << std::endl;
	auto parsedCircuitInfo = readCircuit(circ);
	auto [numQudits, dims, original_circuit] = parsedCircuitInfo;

    Circuit noisyCircuit = original_circuit;

	std::cout << "simsimsimsimsims"<< std::endl;
    if (py::none() == noiseModel) {
        py::dict noiseModelDict = noiseModel.attr("quantum_errors").cast<py::dict>();
        NoiseModel newNoiseModel = parse_noise_model(noiseModelDict);
        printNoiseModel(newNoiseModel);

        std::cout << "======================================================"<< std::endl;

        printCircuit(std::get<2>(parsedCircuitInfo));
        Circuit noisyCircuit = generateCircuit(parsedCircuitInfo, newNoiseModel);

        std::cout << "===================NOISEEEEEEEEE======================="<< std::endl;

        printCircuit(noisyCircuit);
    }


	std::cout << "===================DD SIMULATION======================"<< std::endl;

    //std::tuple<std::basic_string<char>, bool, std::vector<int>, std::basic_string<char>,
    // std::vector<int>, pybind11::object, std::tuple<std::vector<dd::QuantumRegister>, std::vector<dd::Control::Type>>>
    //tag, dag, dims, gate_type, target_qudits, params, control_set
    /*
    py::list paramsNew;

    // Convert field_0 and field_1 to integers (assuming they are integers)
    pybind11::int_ leva = py::cast<int>(0);
    pybind11::int_ levb = py::cast<int>(1);
    pybind11::int_ c = py::cast<int>(1);
    paramsNew.append(leva);
    paramsNew.append(levb);
    paramsNew.append(c);
    paramsNew.append(py::float_(0.));

    // Construct the tuple using the "Instruction" alias
    Instruction first = std::make_tuple(
            "csum",                    // string
            false,                  // bool
            std::vector<int>{3,2},    // vector<int>
            "TWO",               // string
            std::vector<int>{1,0},    // vector<int>
            paramsNew,              // py::object
            std::make_tuple(
                    std::vector<dd::QuantumRegister>{},  // vector<dd::QuantumRegister>
                    std::vector<dd::Control::Type>{}     // vector<dd::Control::Type>
                    )
    );
    Circuit circuit = {first};*/
    std::cout << "op get in "<< std::endl;
    CVec myList = ddsimulator(static_cast<dd::QuantumRegisterCount>(numQudits), static_cast<std::vector<size_t>>(dims), noisyCircuit);

    std::cout << "state vec 0  "<< std::endl;
    for (const auto& elem : myList) {
        std::cout << elem << std::endl;
    }

    py::list result = complex_vector_to_list(myList);

    std::cout << "state vec 1  "<< std::endl;

    return result;

}

PYBIND11_MODULE(pymisim, m) {

    m.doc() = "Example module"; // Optional module docstring
    m.def("state_vector_simulation", &stateVectorSimulation, "state_vector_simulation");
}

