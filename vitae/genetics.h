#pragma once

#include <vector>
#include <tuple>

#include <algorithm>
#include <numeric>
#include <random>
#include <assert.h>

//#define LOG_CYCLE

template<int... Is>
struct seq { };

template<int N, int... Is>
struct gen_seq : gen_seq<N - 1, N - 1, Is...> { };

template<int... Is>
struct gen_seq<0, Is...> : seq<Is...> { };

struct environment_params {
	float desired_fitness;
	float reproduction_pairs_ratio;
	float successful_mate_threshold; //How 'picky' the genes are in selecting a partner <0..1> 0 - chooses any, 1 - only the same or better. Of course you can set it to value greater than 1. And then the real 'pickiness' begins...
	size_t max_generation_count;
	size_t max_population;
	size_t max_gene_age;
	size_t max_num_of_children_per_pair;
	size_t mutation_chance; // should be between 0 and 100
	bool enable_aging;
};

template <typename T>
class gene {
public:
	gene(std::vector<T> input, float fit, size_t gen_size, T min_v, T max_v) : data(input) { fitness = fit; gene_size = gen_size; min_val = min_v; max_val = max_v; }
	float get_fitness() const { return fitness; }
	const std::vector<T>& get_data() const { return data;  }
	void set_fitness(float new_fit) { fitness = new_fit; }
	void inc_age() { ++age; }
	size_t get_size() const { return gene_size; }
	size_t get_age() const { return age; } 
	bool modify_gene_field(size_t index, T value) {
		if (index >= data.size())
			return false;
		data[index] = value;
		return true;
	}
	int get_min_bound() const { return min_val; }
	int get_max_bound() const { return max_val; }
private:
	std::vector<T> data;
	float fitness;
	size_t gene_size;
	size_t age = 0;
	int min_val;
	int max_val;
};

template <typename T, class... Types>
struct fitness_function {
	fitness_function(float(*n_f)(const gene<T>&, std::tuple<Types...>), std::tuple<Types...> params) { f = n_f; args = params; }
	float(*f)(const gene<T>&, std::tuple<Types...>) = nullptr;
	std::tuple<Types...> args;
};


template <typename T, class... Types>
class genes_environment {
public:
	genes_environment(std::vector<gene<T>> entry_population, environment_params& parameters, fitness_function<T, Types...> func_params) : population(entry_population), params(parameters), func(func_params) {
		for (auto& current_gene : population)
			current_gene.set_fitness(func.f(const_cast<const gene<T>&>(current_gene), func.args ) );
		std::sort(population.begin(), population.end(), [](const auto& a, const auto&b) {return a.get_fitness() > b.get_fitness(); });
	}
	bool step_simulation(size_t seed) {
		if (func.f == nullptr)
			return false;
		if (population.empty())
			return false;
		if (generation >= params.max_generation_count)
			return false;
		if (population[0].get_fitness() >= params.desired_fitness )
			return false;

		const auto& max_age = params.max_gene_age;
		const auto& desired_fitness = params.desired_fitness;
		
		std::default_random_engine rand(seed);

		for (auto& gene : population) {
			//incerement age
			gene.inc_age();
		}
		//remove gene if it's too old (weak genes age faster)

#ifdef  LOG_CYCLE
		size_t old_pop = population.size();
#endif //  LOG_CYCLE

		if (params.enable_aging)
			population.erase(std::remove_if(population.begin(), population.end(), [&max_age, &desired_fitness](auto& gene){
				return gene.get_age() * desired_fitness / gene.get_fitness() >= max_age;
			}), population.end());

#ifdef  LOG_CYCLE
		std::cout << old_pop - population.size() << " have passed out!\n";
		old_pop = population.size();
#endif //  LOG_CYCLE

		//if the population excesses limit remove the weakest genes
		if(population.size() > params.max_population)
			population.erase(population.end() - population.size() + params.max_population + 1, population.end() );
		float total_fitness = std::accumulate(population.begin(), population.end(), 0.0f, [](float& a, const gene<T>& b) {return a + b.get_fitness();});
#ifdef  LOG_CYCLE
		std::cout << old_pop - population.size() << " have died due to overpopulation\nTotal population fitness: " << total_fitness << "\n" << population.size() * params.reproduction_pairs_ratio << " pairs have to be bounded.\n";
#endif //  LOG_CYCLE
		//choose genes to reproduce themselves(the youngests have priority)
		size_t assigned_pairs = 0;
		std::vector<gene<T>> bred_genes;
		while (assigned_pairs <= population.size() * params.reproduction_pairs_ratio) {
			for (auto it = population.begin(); it != population.end(); ++it) {
				if (rand() % static_cast<int>(ceil(total_fitness) ) <= ceil(it->get_fitness())) {
					const auto& x = *it; //const reference to the gene that succeed to reproduce.
					size_t chosen_partner_id = 0;
					while (true) {
						//random the second parner untill it fulfills the gene's criteria
						chosen_partner_id = rand() % population.size(); //There's a 1/population.size() chance that the gene chooses itself. Why should we to iterrupt it, anyway?
						if (population[chosen_partner_id].get_fitness() >= x.get_fitness() * params.successful_mate_threshold)
							break;
					}
					const auto& y = population[chosen_partner_id];
					size_t offsprint_count = rand() % params.max_num_of_children_per_pair + 1;
#ifdef  LOG_CYCLE
					std::cout << it - population.begin() << " (" << it->get_fitness() << ") have chosen " << chosen_partner_id << " (" << population[chosen_partner_id].get_fitness() << ") and are about to have " << offsprint_count << " children\n"; 
#endif //  LOG_CYCLE
					for(size_t i = 0; i < offsprint_count; ++i)
						bred_genes.push_back(reproduce(x, y, rand));
					++assigned_pairs;
				}
			}
		}

#ifdef  LOG_CYCLE
		std::cout << bred_genes.size() << " new genes were born.\n";
#endif //  LOG_CYCLE
		population.insert(population.end(), bred_genes.begin(), bred_genes.end());
		std::sort(population.begin(), population.end(), [](const auto& a, const auto&b) {return a.get_fitness() > b.get_fitness(); });
		++generation;
		return true;
	};
	const std::vector<gene<T> >& get_population() { return population; }
	size_t get_current_generation() { return generation; }

	gene<T> reproduce(const gene<T> &x, const gene<T> &y, std::default_random_engine& rand) {
		assert(x.get_data().size() == y.get_data().size());
		size_t cutting_point = rand() % x.get_size(); 
		std::vector<T> data(x.get_data().size());
		std::copy(x.get_data().begin(), x.get_data().begin() + cutting_point + 1, data.begin());
		std::copy(y.get_data().begin() + cutting_point + 1, y.get_data().end(), data.begin() + cutting_point + 1);
		if (rand() % 100 <= params.mutation_chance)
			data[rand() % data.size()] = rand() % (x.get_max_bound() - x.get_min_bound()) + x.get_min_bound();

		auto new_gene = gene<T>(data, 0.0f, x.get_size(), x.get_min_bound(), x.get_max_bound() );
		new_gene.set_fitness(func.f(const_cast<const gene<T>&>(new_gene), func.args) );
		return new_gene;
	}

private:
	std::vector<gene<T> > population;
	size_t generation = 0;
	environment_params& params;
	fitness_function<T, Types...> func;

};

template <typename T>
std::vector<gene<T>> random_entry_population(size_t count, size_t genes_size, int lower_bound, int higher_bound, size_t seed ) {
	assert(higher_bound > lower_bound); //Hola, hola!
	std::vector<gene<T>> output;
	std::default_random_engine rand(seed);
	for (size_t i = 0; i < count; ++i) {
		std::vector<T> current;
		for (size_t j = 0; j < genes_size; ++j)
			current.push_back(rand() % (higher_bound - lower_bound) + lower_bound);
		output.push_back(gene<T>(current, 0.0f, genes_size, lower_bound, higher_bound) );
	}
	return output;
}