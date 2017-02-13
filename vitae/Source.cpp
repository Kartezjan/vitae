#include <iostream>
#include <string>
#include <ctime>

#include "genetics.h"

float horizontal_search(const std::vector<bool>& chessboard, const size_t& size, const size_t& x) {
	float result = 0.0f;
	for (size_t i = 0; i < size; ++i) {
		if (chessboard[x*size + i])
			result += 1.0f;
	}
	return result - 1.0f;
}

float vertical_search(const std::vector<bool>& chessboard, const size_t& size, const size_t& y) {
	float result = 0.0f;
	for (size_t i = 0; i < size; ++i) {
		if (chessboard[i*size + y])
			result += 1.0f;
	}
	return result - 1.0f;
}

float diagonal_search(const std::vector<bool>& chessboard, const size_t& size, const int& x, const int& y) {
	float result = 0.0f;
	int i = -1, j = 1;
	while (x + i >= 0 && y + j < size) {
		if (chessboard[(x + i) * size + (y + j)])
			result += 1.0f;
		--i, ++j;
	}

	i = 1, j = 1;
	while (x + i < size && y + j < size) {
		if (chessboard[(x + i) * size + (y + j)])
			result += 1.0f;
		++i, ++j;
	}

	i = 1, j = -1;
	while (x + i < size && y + j >= 0) {
		if (chessboard[(x + i) * size + (y + j)])
			result += 1.0f;
		++i, --j;
	}

	i = -1, j = -1;
	while (x + i >= 0 && y + j >= 0) {
		if (chessboard[(x + i) * size + (y + j)])
			result += 1.0f;
		--i, --j;
	}

	return result;
}

float calc_fitness(const std::vector<bool>& chessboard, const size_t& size, const std::vector<char>& queen_positions) {
	float fitness = 60.f;
	for (size_t i = 0; i < size; ++i) {
		fitness -= horizontal_search(chessboard, size, i);
		fitness -= vertical_search(chessboard, size, queen_positions.data()[i]);
		fitness -= diagonal_search(chessboard, size, i, queen_positions.data()[i]);
	}
	return fitness;
}

void print_solution(const std::vector<char>& solution, size_t row_len) {
	std::vector<bool> chessboard(row_len*row_len);
	std::fill(chessboard.begin(), chessboard.end(), 0);
	for (size_t i = 0; i < row_len; ++i) {
		chessboard[i * 8 + solution[i]] = true;
	}
	for (size_t i = 0; i < chessboard.size(); ++i) {
		if (i % row_len == 0)
			std::cout << "\n" << std::string(4 * row_len, '-') << "\n";
		if (chessboard[i])
			std::cout << "|X| ";
		else
			std::cout << "| | ";
	}
	std::cout << "\n" << std::string(4 * row_len, '-') << "\n";
}

int main(void) {

	size_t size = 8;
	std::vector<bool> chessboard(size*size);

	environment_params params{ 60.0f, 0.8f, 0.8f, 100, 10000, 20, 2, 5, false };
	fitness_function<char, std::vector<bool>*, size_t* > func_params([](const gene<char> &copy, std::tuple<std::vector<bool>*, size_t*> tuple) -> float {
		const size_t& size = *std::get<size_t*>(tuple);
		std::vector<bool>& chessboard = *std::get<std::vector<bool>*>(tuple);
		std::fill(chessboard.begin(), chessboard.end(), 0);

		for (size_t i = 0; i < size; ++i) {
			chessboard[i * 8 + copy.get_data()[i]] = true;
		}
		return calc_fitness(chessboard, size, copy.get_data() );
	}, std::make_tuple<std::vector<bool>*, size_t* >(&chessboard, &size));
	genes_environment<char, std::vector<bool>*, size_t*> engine(random_entry_population<char>(100, 8, 0, 8, 0x1337), params, func_params);

	do {
		std::cout << "current generation: " << engine.get_current_generation() << "\n" <<
		std::string(30, '-') << "\n";
		if (!engine.get_population().empty())
			std::cout << "population: " << engine.get_population().size() << "\n" <<
			"max fitness: " << std::max_element(engine.get_population().begin(), engine.get_population().end(), [](const auto& a, const auto& b) {return a.get_fitness() < b.get_fitness();})->get_fitness() << "\n" <<
			"min fitnesss: " << std::min_element(engine.get_population().begin(), engine.get_population().end(), [](const auto& a, const auto& b) {return a.get_fitness() < b.get_fitness();})->get_fitness() << "\n";
		else
			std::cout << "All genes are dead!\n";
		std::cout << "\n";
	} while (engine.step_simulation(time(NULL)));

	std::cout << "\n\n\nBest solution:\n";
	print_solution(engine.get_population()[0].get_data(), size);
	return 0;
}