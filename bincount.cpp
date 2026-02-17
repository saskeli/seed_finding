/*
 * Copyright (c) 2026 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cstddef>
#include <iostream>
#include <string>
#include <unordered_map>


int main(int argc, char **argv)
{
	std::ios_base::sync_with_stdio(false);
	
	std::unordered_map <std::string, std::size_t> counts;
	std::string buffer;
	while (std::getline(std::cin, buffer))
		++counts[buffer];
	
	for (auto const &kv : counts)
		std::cout << kv.second << '\t' << kv.first << '\n';
	
	return 0;
}
