#pragma once

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <memory>

#include "ATKConnectorDll.h"

const std::string space = " ";
const double TwoDays = 86400.0 * 2.0;

/*------------------------------------------------------------------------------*
 * Data of one satellite :														*
 * 1. Intiial OE : a0, e0, i0, Omega0, omega0, f0								*
 * 2. Impulse time list															*
 * 3. Impulse list : dvx, dvy, dvz												*
 *------------------------------------------------------------------------------*/
class SatelliteData {

private:
	std::vector<double> initial_coe_;																		// Initial orbit elements (km, s)
	std::vector<double> impulse_t_list_;																	// Impulse time list (s)
	std::vector<std::vector<double>> impulse_list_;														// Impulse list (km, s)
	
// Constructor and destructor
public:
	SatelliteData() = default;
	~SatelliteData() = default;

// Functions: set value and get value, in order to protect the data
public:
	// Set value
	void set_value(const std::vector<double>, const std::vector<double>, const std::vector<std::vector<double>>);

	// Get value
	std::vector<double> get_initial_coe() { return this->initial_coe_; }
	std::vector<double> get_impulse_t_list() { return this->impulse_t_list_; }
	std::vector<std::vector<double>> get_impulse_list() { return this->impulse_list_; }
	
};


/*------------------------------------------------------------------------------*
 * Data of the final result :													*
 * 1. Read the data and Push back the satellites								*
 * 2. Write the information into ATK											*
 * TXT format:																	*
 * 1 a0 e0 i0 Omega0 omega0 f0													*
 * 1 t1 dv1x dv1y dv1z															*
 * 2 a0 e0 i0 Omega0 omega0 f0													*
 * 2 t1 dv1x dv1y dv1z															*
 *------------------------------------------------------------------------------*/
class ResultData {

// Satellites
private:
	std::vector<SatelliteData*> satellites_;
	int conID_;

// Constructor and destructor
public:
	ResultData() = default;
	~ResultData() { std::vector<SatelliteData*>().swap(satellites_); }

// Read txt and push back satellites
// Write into ATK
public:
	void read_data(const std::string);
	void write_atk();

private:
	void add_sat(SatelliteData* sat) { this->satellites_.push_back(sat); }
	void write_single_sat(SatelliteData*, const int, const int);
	void create_sateliite(const std::string, std::string&) const;
	void add_initialstate(const std::vector<double>, const std::string) const;
	void add_propagator(const double, const std::string, const std::string) const;
	void add_impulse(const std::vector<double>, const std::string, const std::string) const;
};




