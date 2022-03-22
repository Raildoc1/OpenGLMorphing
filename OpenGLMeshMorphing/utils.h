#pragma once
#include <ctime>
#include <string>
#include <iostream>
#include <filesystem>

void print_time_stamp(const clock_t& prev, const std::string msg) {
	std::cout << std::setw(30) << msg << std::setw(30) << float(clock() - prev) / CLOCKS_PER_SEC << " seconds." << std::endl;
}