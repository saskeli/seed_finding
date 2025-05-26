# Copyright (c) 2025 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).

import argparse


def defined_characters(ss: str) -> int:
	retval = 0
	for cc in ss:
		if '.' != cc:
			retval += 1
	return retval


def cmp(c1: str, c2: str) -> bool:
	if '.' == c1:
		return False
	if '.' == c2:
		return False
	return c1 == c2


def alignment_score(s1: str, s2: str, offset: int) -> int:
	retval = 0
	if offset < 0:
		assert -offset < len(s1)
		limit = min(len(s1) + offset, len(s2))
		for ii in range(0, limit):
			retval += cmp(s1[ii - offset], s2[ii])
	else:
		assert offset < len(s1)
		limit = min(len(s2) - offset, len(s1))
		for ii in range(0, limit):
			retval += cmp(s1[ii], s2[ii + offset])
	return retval


def max_alignment_score(s1: str, s2: str) -> tuple[int, int]:
	score = 0
	offset = 0
	ll = len(s1)
	if 0 == ll:
		return 0, 0

	for ii in range(-ll + 1, ll):
		current_score = alignment_score(s1, s2, ii)
		if score < current_score:
			score = current_score
			offset = ii

	return score, offset


def huddinge_distance(s1: str, s2: str) -> tuple[int, int]:
	dd = max(defined_characters(s1), defined_characters(s2))
	score, offset = max_alignment_score(s1, s2)
	return dd - score, offset


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = "Calculate the Huddinge distance of the given strings.")
	parser.add_argument('s1', type = str, help = "First string")
	parser.add_argument('s2', type = str, help = "Second string")
	args = parser.parse_args()

	distance, offset = huddinge_distance(args.s1, args.s2)
	print(f"Distance: {distance} offset: {offset}")
