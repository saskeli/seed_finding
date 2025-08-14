from abc import ABC, abstractmethod
import argparse
import concurrent.futures
from dataclasses import dataclass, field
import datetime
import edlib
import gzip
import math
import parse
import os
import sys
import typing


INPUT_PREFIX: typing.Final[str] = "/home/tnorri/selex/seed-finder-output-2"


def make_iupac_equality_list():
	# IUPAC characters are allowed in the expected seed.
	equalities = [
		# Expected on the left.
		("R", "AG"),
		("Y", "CT"),
		("S", "GC"),
		("W", "AT"),
		("K", "GT"),
		("M", "AC"),
		("B", "CGT"),
		("D", "AGT"),
		("H", "ACT"),
		("V", "ACG"),
		("N", "ACGT"), # “.” on the following line.
		("ACGTRYSWKMBDHVN", ".")
	]

	def helper():
		for lhss, rhss in equalities:
			for cc in lhss:
				for cc_ in rhss:
					yield cc_, cc # Expected seed on the right.
	return list(helper())

IUPAC_EQUALITIES: typing.Final[[typing.Tuple[str, str]]] = make_iupac_equality_list()


COUNT_PARSER: typing.Final[parse.Parser] = parse.compile(r'({:d}, {:d})')


def reverse_complement_expected_seed(expected_seed):
	"Build the reverse complement of the given expected seed"

	equalities = {
		'A': 'T',
		'C': 'G',
		'G': 'C',
		'T': 'A',
		'R': 'Y',
		'Y': 'R',
		'S': 'S',
		'W': 'W',
		'K': 'M',
		'M': 'K',
		'B': 'V',
		'D': 'H',
		'H': 'D',
		'V': 'B',
		'N': 'N'
	}

	def helper():
		for cc in reversed(expected_seed):
			yield equalities[cc]

	return ''.join(helper())


class Warning(Exception):
	pass


class Error(Exception):
	pass


def log(message: str):
	current_time = datetime.datetime.now()
	current_time_ = current_time.strftime("%H:%M:%S")
	print(f"[{current_time_}] {message}", file = sys.stderr)


def log_warning(message: str):
	log(f"WARNING: {message}")


def log_error(message: str):
	log(f"ERROR: {message}")


def normalised_edit_similarity(distance, mm):
	# From Lopresti, Zhou: Retrieval Strategies for Noisy Text, eq. 6.
	return 1.0 / math.exp(distance / (mm - distance))


@dataclass
class Task(ABC):
	identifier: str
	expected_seed: str
	should_test_reverse_complement: bool
	input_limit: int
	should_read_tail: bool

	@abstractmethod
	def align(self, *, seed: str, expected_seed: str, is_reverse_complement: bool, counts: typing.Tuple[int, int], p_val: float, priority: float):
		pass

	@abstractmethod
	def output(self):
		pass

	def read_from_file(self, path: str):
		with open_(path) as fp:
			for line in fp:
				if line.startswith("#"):
					continue
				yield line

	def read_head(self, path: str):
		for lineno, line in enumerate(self.read_from_file(path)):
			if self.input_limit <= lineno:
				return
			yield line

	def read_tail(self, path: str):
		# Create a circular buffer for the lines.
		lines = [""] * self.input_limit
		count = 0
		current_idx = 0
		for lineno, line in enumerate(self.read_from_file(path), start = 1):
			lines[current_idx] = line
			current_idx = lineno % self.input_limit
			count = lineno
		if self.input_limit < count:
			yield from iter(lines[current_idx:])
		yield from iter(lines[:current_idx])

	def input_lines(self, path: str):
		if -1 == self.input_limit:
			yield from self.read_from_file(path)
		else:
			if self.should_read_tail:
				yield from self.read_tail(path)
			else:
				yield from self.read_head(path)

	def run(self):
		seed_finder_output_path = f"{INPUT_PREFIX}/{self.identifier}.tsv.gz"

		try:
			expected_seed_rc = reverse_complement_expected_seed(self.expected_seed) if self.should_test_reverse_complement else ""
			for line_ in self.input_lines(seed_finder_output_path):
				line = line_.rstrip("\n")
				try:
					seed, counts_, p_val, priority = line.split("\t")
					res = COUNT_PARSER.parse(counts_)
					if res is None:
						raise Error(f'Unexpected line format: "{line}"')
					assert isinstance(res, parse.Result)
					counts = res.fixed
				except ValueError:
					raise Error(f'Unexpected line format: "{line}"')

				def align_(seed: str, expected_seed: str, is_reverse_complement: bool):
					self.align(
						seed = seed,
						expected_seed = expected_seed,
						is_reverse_complement = is_reverse_complement,
						counts = counts,
						p_val = float(p_val),
						priority = float(priority)
					)

				align_(seed, self.expected_seed, False)
				if self.should_test_reverse_complement:
					align_(seed, expected_seed_rc, True)
		except OSError as exc:
			raise Warning(f"Skipping {seed_finder_output_path}: {exc}")


@dataclass
class PathAlignmentTask(Task):
	results: typing.List[typing.Tuple[bool, str, str, str, str]] = field(default_factory = list)

	def align(self, *, seed: str, expected_seed: str, is_reverse_complement: bool, counts: typing.Tuple[int, int], p_val: float, priority: float):
		res = edlib.align(seed, expected_seed, additionalEqualities = IUPAC_EQUALITIES, mode = "NW", task = "path")
		aa = edlib.getNiceAlignment(res, seed, expected_seed)
		self.results.append((is_reverse_complement, seed, aa["target_aligned"], aa["matched_aligned"], aa["query_aligned"]))

	def output(self):
		for (is_reverse_complement, seed, target_aligned, matched_aligned, query_aligned) in self.results:
			expected_seed = reverse_complement_expected_seed(self.expected_seed) if is_reverse_complement else self.expected_seed
			print(f"Identifier:            {self.identifier}")
			print(f"Expected seed:         {expected_seed}")
			print(f"Is reverse complenent: {is_reverse_complement}")
			print(f"Found seed:            {seed}")
			print(target_aligned)
			print(matched_aligned)
			print(query_aligned)


@dataclass
class DistanceAlignmentTask(Task):

	@dataclass
	class Result:
		is_reverse_complement: bool
		seed: str
		edit_distance: int
		normalised_edit_similarity: float
		counts: typing.Tuple[int, int]
		p_val: float
		priority: float

	results: typing.List[Result] = field(default_factory = list)

	def align(self, *, seed: str, expected_seed: str, is_reverse_complement: bool, counts: typing.Tuple[int, int], p_val: float, priority: float):
		res = edlib.align(seed, expected_seed, additionalEqualities = IUPAC_EQUALITIES, mode = "NW", task = "distance")
		ed = res['editDistance']
		ned = normalised_edit_similarity(ed, len(seed))
		self.results.append(self.Result(is_reverse_complement, seed, ed, ned, counts, p_val, priority))

	def output(self):
		def format_bool(value: bool):
			return 1 if value else 0

		for rr in self.results:
			expected_seed = reverse_complement_expected_seed(self.expected_seed) if rr.is_reverse_complement else self.expected_seed
			print(f"{self.identifier}\t{expected_seed}\t{format_bool(rr.is_reverse_complement)}\t{rr.seed}\t{rr.edit_distance}\t{rr.normalised_edit_similarity}\t{rr.counts}\t{rr.p_val}\t{rr.priority}")


def open_(path: str) -> typing.TextIO:
	if ".gz" == path[-3:]:
		return gzip.open(path, mode = "rt")
	else:
		return open(path, "r")


def read_motif_input(fp: typing.TextIO) -> typing.Iterable[typing.Tuple[str, str]]:
	log("Reading motifs from stdin…")
	is_first = True

	for motif_lineno, motif_line_ in enumerate(sys.stdin, start = 1):
		if 0 == motif_lineno % 100:
			log(f"Line {motif_lineno}…")

		motif_line = motif_line_.rstrip("\n")
		fields = motif_line.split("\t")

		if is_first:
			is_first = False

			# Sanity check.
			# FIXME: Since the check is input-dependent, assertions should not be used.
			assert fields[0] == "ID"
			assert fields[7] == "cycle"
			assert fields[8] == "cycle_background"
			assert fields[10] == "seed"
			assert fields[11] == "CSC_SELEX_filename"
			assert fields[12] == "CSC_SELEX_background_filename"
			continue

		identifier = fields[0]
		expected_seed = fields[10]

		yield identifier, expected_seed


def prepare_tasks(it: typing.Iterable[typing.Tuple[str, str]], *, should_output_alignment = False, should_test_reverse_complement = False, input_limit = -1, should_read_tail = False):
	for identifier, expected_seed in it:
		if should_output_alignment:
			task = PathAlignmentTask(identifier, expected_seed, should_test_reverse_complement, input_limit, should_read_tail)
		else:
			task = DistanceAlignmentTask(identifier, expected_seed, should_test_reverse_complement, input_limit, should_read_tail)
		assert task
		yield task


def run_task(task: Task):
	task.run()
	return task


def handle_results(done: typing.Set[concurrent.futures.Future]):
	for ff in done:
		res: typing.Optional[Task] = None
		try:
			res = ff.result()
		except Error as exc:
			log_error(f"{exc}")
			continue
		except Warning as exc:
			log_warning(f"{exc}")
			continue
		except Exception as exc:
			log_error(f"Got an exception: {type(exc)}: {exc}")
			continue
		assert res
		res.output()


def main():
	parser = argparse.ArgumentParser(description = 'Compare seed_finder’s results to expected seeds')
	parser.add_argument("--output-alignment", action = "store_true", help = "Output the alignments instead of alignment scores")
	parser.add_argument("--test-reverse-complement", action = "store_true", help = "Test also the reverse complement of the expected seed")
	parser.add_argument("--input-limit", type = int, default = -1, metavar = "COUNT", help = "Check only the first COUNT found seeds, -1 for no limit")
	parser.add_argument("--tail", action = "store_true", help = "Read the last lines of the input instead of the first ones")
	args = parser.parse_args()

	log(f"Using seed_finder outputs at {INPUT_PREFIX}…")

	if not(args.output_alignment):
		print("#IDENTIFIER\tEXPECTED_SEED\tIS_REVERSE_COMPLEMENT\tSEED\tEDIT_DISTANCE\tNORMALISED_EDIT_SIMILARITY\tCOUNTS\tP_VAL\tPRIORITY")

	max_workers = os.cpu_count() or 1 # FIXME: Change to process_cpu_count().
	# FIXME: I’m not sure what causes the memory usage of the worker processes to increase over time. Apparently there is a leak (or something similar) since limiting the number of tasks per child helps.
	with concurrent.futures.ProcessPoolExecutor(max_workers = max_workers, max_tasks_per_child = 1) as executor:
		max_tasks = 2 * max_workers
		running_tasks = set()
		for tt in prepare_tasks(
			read_motif_input(sys.stdin),
			should_output_alignment = args.output_alignment,
			should_test_reverse_complement = args.test_reverse_complement,
			input_limit = args.input_limit,
			should_read_tail = args.tail
		):
			while max_tasks <= len(running_tasks):
				done, not_done = concurrent.futures.wait(running_tasks, return_when = concurrent.futures.FIRST_COMPLETED)
				handle_results(done)
				running_tasks = not_done
			ff = executor.submit(run_task, tt)
			running_tasks.add(ff)
		done, _ = concurrent.futures.wait(running_tasks, return_when = concurrent.futures.ALL_COMPLETED)
		handle_results(done)


if __name__ == "__main__":
	main()
