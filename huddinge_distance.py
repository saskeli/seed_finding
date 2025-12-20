# Copyright (c) 2025 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).

import argparse


def defined_characters(ss: str, wildcard: str) -> int:
	retval = 0
	for cc in ss:
		if wildcard != cc:
			retval += 1
	return retval


def cmp(c1: str, c2: str, wildcard: str) -> bool:
	if wildcard == c1:
		return False
	if wildcard == c2:
		return False
	return c1 == c2


def alignment_score(s1: str, s2: str, offset: int, wildcard: str) -> int:
	retval = 0
	if offset < 0:
		assert -offset < len(s1)
		limit = min(len(s1) + offset, len(s2))
		for ii in range(0, limit):
			retval += cmp(s1[ii - offset], s2[ii], wildcard)
	else:
		assert offset < len(s1)
		limit = min(len(s2) - offset, len(s1))
		for ii in range(0, limit):
			retval += cmp(s1[ii], s2[ii + offset], wildcard)
	return retval


def max_alignment_score(s1: str, s2: str, wildcard: str) -> tuple[int, int]:
	score = 0
	offset = 0
	ll = len(s1)
	if 0 == ll:
		return 0, 0

	for ii in range(-ll + 1, ll):
		current_score = alignment_score(s1, s2, ii, wildcard)
		if score < current_score:
			score = current_score
			offset = ii

	return score, offset


def huddinge_distance(s1: str, s2: str, wildcard: str) -> tuple[int, int]:
	dd = max(defined_characters(s1, wildcard), defined_characters(s2, wildcard))
	score, offset = max_alignment_score(s1, s2, wildcard)
	return dd - score, offset


def huddinge_distance_at(s1: str, s2: str, offset: int, wildcard: str) -> int:
	dd = max(defined_characters(s1, wildcard), defined_characters(s2, wildcard))
	score = alignment_score(s1, s2, offset, wildcard)
	return dd - score


def len_eq_1(ss: str):
	if 1 == len(ss):
		return ss
	raise argparse.ArgumentTypeError("Value must have length of exactly one")


def main():
	parser = argparse.ArgumentParser(description = "Calculate the Huddinge distance of the given strings.")
	parser.add_argument("--wildcard", type = len_eq_1, default = ".", help = "Wildcard character")
	parser.add_argument("--offset", type = int, default = None, help = "Calculate distance for the given offset")
	parser.add_argument('s1', type = str, help = "First string")
	parser.add_argument('s2', type = str, help = "Second string")
	args = parser.parse_args()

	if args.offset is None:
		distance, offset = huddinge_distance(args.s1, args.s2, args.wildcard)
	else:
		offset = args.offset
		distance = huddinge_distance_at(args.s1, args.s2, offset, args.wildcard)

	print(f"Distance: {distance} offset: {offset}")
	print("Alignment:")
	offset_ = abs(offset)
	spaces = abs(offset_) * " "
	if 0 <= offset:
		print(f"{spaces}{args.s1}")
		print(args.s2)
	else:
		print(args.s1)
		print(f"{spaces}{args.s2}")
	alignment = "".join(
		map(
			lambda x: "*" if cmp(x[0], x[1], args.wildcard) else " ",
			zip(args.s1[(0 if 0 <= offset else offset_):], args.s2[(0 if offset < 0 else offset_):])
		)
	)
	print(f"{spaces}{alignment}")


if __name__ == "__main__":
	main()
