#include "ResultData.h"



// If elements exist, refuse to modify the data
void SatelliteData::set_value(const std::vector<double> coe, const std::vector<double> t_list, const std::vector<std::vector<double>> dv_list) {
	if (!this->initial_coe_.empty()) {
		std::cout << "Initial orbit elements exist, refuse to set new value." << std::endl;
	} else if (!this->impulse_t_list_.empty()) {
		std::cout << "Impulse moments list exists, refuse to set new value." << std::endl;
	} else if (!this->impulse_list_.empty()) {
		std::cout << "Impulse list exists, refuse to set new value." << std::endl;
	} else {
		for (int i = 0; i < coe.size(); i++) this->initial_coe_.push_back(coe[i]);
		for (int i = 0; i < t_list.size(); i++) this->impulse_t_list_.push_back(t_list[i]);
		for (int i = 0; i < dv_list.size(); i++) this->impulse_list_.push_back(dv_list[i]);
	}
}


void ResultData::read_data(const std::string filename) {
	std::ifstream file(filename);
	std::string line;

	while (std::getline(file, line)) {
		std::istringstream iss(line);
		std::vector<double> initial_coe(6);
		for (double& val : initial_coe) iss >> val;

		std::vector<double> impulse_t_list;
		std::vector<std::vector<double>> impulse_list;
		while (std::getline(file, line)) {
			std::istringstream iss_impulse(line);
			std::vector<double> impulse_data;
			double value;
			while (iss_impulse >> value) impulse_data.push_back(value);

			if (impulse_data.size() == 4) {
				impulse_t_list.push_back(impulse_data[0]);
				impulse_list.push_back({ impulse_data[1], impulse_data[2], impulse_data[3] });
			} else if (impulse_data.size() == 6) {
				// Push the current satellite data to the vector
				SatelliteData* sat = new SatelliteData();
				sat->set_value(initial_coe, impulse_t_list, impulse_list);
				add_sat(sat);

				// Start reading the next satellite data
				initial_coe = impulse_data;
				impulse_t_list.clear();
				impulse_list.clear();
			}
		}

		SatelliteData* sat = new SatelliteData();
		sat->set_value(initial_coe, impulse_t_list, impulse_list);
		add_sat(sat);
	}
}

void ResultData::write_atk() {
	conID_ = atkOpen();
	for (auto i = 0; i < satellites_.size(); i++)
		this->write_single_sat(satellites_[i], conID_, i+1);
	atkClose(conID_);
}


// If the satellite does not exist, then create one
void ResultData::create_sateliite(const std::string id, std::string& sat_path) const {
	atkConnect(conID_, "New", "/ Satellite Satellite" + id);														// New satellite
	atkConnect(conID_, "New", "/ */Satellite/Satellite" + id + "/Sensor Sensor" + id);							// New sensor
	atkConnect(conID_, "Define", "*/Satellite/Satellite" + id + "/Sensor/Sensor" + id + "Rectangular 20 20");	// Sensor parameter
	sat_path =  "*/Satellite/Satellite" + id;
}


// Add initial state of the satellite
void ResultData::add_initialstate(const std::vector<double> initial_coe, const std::string sat_path) const{
	// Set ATK commands
	std::string set_initialstate_command = sat_path + space + "SetValue MainSequence.SegmentList.Initial_State.CoordinateType Keplerian";
	std::string add_initialstate_sma_command = sat_path + space + "SetValue MainSequence.SegmentList.Initial_State.InitialState.Keplerian.sma" + space + std::to_string(initial_coe[0]) + space + "km";
	std::string add_initialstate_ecc_command = sat_path + space + "SetValue MainSequence.SegmentList.Initial_State.InitialState.Keplerian.ecc" + space + std::to_string(initial_coe[1]);
	std::string add_initialstate_inc_command = sat_path + space + "SetValue MainSequence.SegmentList.Initial_State.InitialState.Keplerian.inc" + space + std::to_string(initial_coe[2]) + space + "rad";
	std::string add_initialstate_RAAN_command = sat_path + space + "SetValue MainSequence.SegmentList.Initial_State.InitialState.Keplerian.RAAN" + space + std::to_string(initial_coe[3]) + space + "rad";
	std::string add_initialstate_w_command = sat_path + space + "SetValue MainSequence.SegmentList.Initial_State.InitialState.Keplerian.w" + space + std::to_string(initial_coe[4]) + space + "rad";
	std::string add_initialstate_ta_command = sat_path + space + "SetValue MainSequence.SegmentList.Initial_State.InitialState.Keplerian.ta" + space + std::to_string(initial_coe[5]) + space + "rad";

	// Add the data
	atkConnect(conID_, "Astrogator", set_initialstate_command);
	atkConnect(conID_, "Astrogator", add_initialstate_sma_command);
	atkConnect(conID_, "Astrogator", add_initialstate_ecc_command);
	atkConnect(conID_, "Astrogator", add_initialstate_inc_command);
	atkConnect(conID_, "Astrogator", add_initialstate_RAAN_command);
	atkConnect(conID_, "Astrogator", add_initialstate_w_command);
	atkConnect(conID_, "Astrogator", add_initialstate_ta_command);
}

// Add a new propagator
void ResultData::add_propagator(const double duration_time, const std::string sat_path, const std::string prop_id) const {
	atkConnect(conID_, "Astrogator", sat_path + space + "InsertSegment MainSequence.SegmentList.- Propagate");
	atkConnect(conID_, "Astrogator", sat_path + space + "SetValue MainSequence.SegmentList.Propagate" + prop_id + ".StoppingConditions Duration");
	atkConnect(conID_, "Astrogator", sat_path + space + "SetValue MainSequence.SegmentList.Propagate" + prop_id + ".StoppingConditions.Duration.TripValue" + space + std::to_string(duration_time) + space + "sec");
}


// Add impulses
void ResultData::add_impulse(const std::vector<double> impulse_cartesian, const std::string sat_path, const std::string manuv_id) const {
	double dvx = impulse_cartesian[0];
	double dvy = impulse_cartesian[1];
	double dvz = impulse_cartesian[2];

	// Add impulse, if this is the first impulse, there is no maneuver id, or else there should be one
	atkConnect(conID_, "Astrogator", sat_path + space + "InsertSegment MainSequence.SegmentList.- Maneuver");

	// Set the commands
	std::string set_impulseAxis = sat_path + space + "SetValue MainSequence.SegmentList.Maneuver" + manuv_id + ".ImpulsiveMnvr.ThrustAxes Satellite_VNC(Earth)";
	std::string set_impulseValueX = sat_path + space + "SetValue MainSequence.SegmentList.Maneuver" + manuv_id + ".ImpulsiveMnvr.Cartesian.X" + space + std::to_string(dvx) + "km/sec";
	std::string set_impulseValueY = sat_path + space + "SetValue MainSequence.SegmentList.Maneuver" + manuv_id + ".ImpulsiveMnvr.Cartesian.Y" + space + std::to_string(dvy) + "km/sec";
	std::string set_impulseValueZ = sat_path + space + "SetValue MainSequence.SegmentList.Maneuver" + manuv_id + ".ImpulsiveMnvr.Cartesian.Z" + space + std::to_string(dvz) + "km/sec";

	// Set the coordinates of maneuver
	atkConnect(conID_, "Astrogator", set_impulseAxis);
	// X,Y,Z components
	atkConnect(conID_, "Astrogator", set_impulseValueX);
	atkConnect(conID_, "Astrogator", set_impulseValueY);
	atkConnect(conID_, "Astrogator", set_impulseValueZ);
}

// Wrtie the data of one satellite into ATK
void ResultData::write_single_sat(SatelliteData* sat, const int conID, const int sat_id) {
	std::string id = std::to_string(sat_id);
	std::string sat_path;

	//Get data of the satellite
	std::vector<double> impulse_t_list = sat->get_impulse_t_list();
	std::vector<std::vector<double>> impulse_list = sat->get_impulse_list();
	std::vector<double> initial_coe = sat->get_initial_coe();

	// Create a new satellite
	create_sateliite(id, sat_path);
	std::cout << id << space << sat_path << std::endl;

	// Set the mode into Setprop
	atkConnect(conID, "Astrogator", sat_path +  space + "SetProp");

	// Set the initial state
	add_initialstate(initial_coe, sat_path);

	// First propagating segment
	double duration_time = impulse_t_list[0];
	add_propagator(duration_time, sat_path, "");

	for (int i = 0; i < impulse_t_list.size(); i++) {
		if (i != 0) {
			std::string manuv_id = std::to_string(i);
			add_impulse(impulse_list[i], sat_path, manuv_id);
		} else {
			std::string manuv_id = "";
			add_impulse(impulse_list[i], sat_path, manuv_id);
		}

		// Add propagator, if this is the last impulse, the duration time should be t_final - t_imp
		if (i != impulse_t_list.size()-1) {
			std::string prop_id = std::to_string(i+1);
			duration_time = impulse_t_list[i+1] - impulse_t_list[i];
			add_propagator(duration_time, sat_path, prop_id);
		}
		else {
			std::string prop_id = std::to_string(i+1);
			duration_time = TwoDays - impulse_t_list[i];
			add_propagator(duration_time, sat_path, prop_id);
		}
	}

	atkConnect(conID_, "Astrogator", sat_path + space + "RunMCS");
}
