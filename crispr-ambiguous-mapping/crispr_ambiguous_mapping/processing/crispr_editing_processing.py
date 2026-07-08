from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Optional, DefaultDict, Union, Tuple, List, Any
import pandas as pd
from typing import Counter as CounterType
from ..models.mapping_models import GeneralMappingInferenceDict, GeneralAlleleDict, GeneralAlleleCountSeriesDict
from ..models.mapping_models import CompleteInferenceMatchResult, MatchSetSingleInferenceMatchResult, InferenceResult
from ..models.editing_models import (MatchSetWhitelistReporterObservedSequenceCounterSeriesResults, 
                                     MatchSetWhitelistReporterObservedSequenceMutationProfiles, 
                                     LinkedMutationCounters, 
                                     ObservedSequenceMutationProfile)

# Function that has conditional on whether a match result is not an errror
def check_match_result_non_error(match_result):
    return False if match_result is None else match_result.error is None # If match_result is None, treat as error. If match_result is not None, but error is None, then non_error

# Filter dict with observed sequence inference results for only those that do not contain any mapping errors
def get_non_error_dict(observed_guide_reporter_umi_counts_inferred: Union[GeneralMappingInferenceDict, DefaultDict[str, GeneralMappingInferenceDict]], attribute_name: str) -> MatchSetWhitelistReporterObservedSequenceCounterSeriesResults:
    """
    Return a dict of all entries that are non-error.
    Works with both sample-barcode (nested dict) and no-sample-barcode (flat dict) cases.
    """
    non_error_dict = {}

    for key, val in observed_guide_reporter_umi_counts_inferred.items():
        if hasattr(val, "inferred_value"):  # flat dict
            attr = getattr(val.inferred_value, attribute_name, None)
            if check_match_result_non_error(attr):
                non_error_dict[key] = val
        else:  # nested dict (sample-barcode case)
            for inner_key, inner_val in val.items():
                attr = getattr(inner_val.inferred_value, attribute_name, None)
                if check_match_result_non_error(attr):
                    non_error_dict[(key, inner_key)] = inner_val  # tuple key: (sample, guide)

    return non_error_dict

#
# Given the datastructure containing the inference results "observed_guide_reporter_umi_counts_inferred", iterate through the entire datastructure to generate
# another datastructure that contains the observed alleles (protospacer/surrogate/barcode) for each whitelist reporter in either a dictionary or dataframe format.
# This is the foundational datastructure used for analyzing the mutations for each whitelist reporter.
#
def _require_inference_dict(observed_guide_reporter_umi_counts_inferred, caller: str):
    """Raise a clear error if the slim (retain_inference_results=False)
    result's inference dict is None. Saves users from a cryptic AttributeError
    deep inside the allele-building loop."""
    if observed_guide_reporter_umi_counts_inferred is None:
        raise ValueError(
            f"{caller} requires `observed_guide_reporter_umi_counts_inferred` but the "
            f"result object was built with `retain_inference_results=False` (the default). "
            f"Re-run the mapping call with `retain_inference_results=True` to enable "
            f"allele / mutation post-processing."
        )


def get_matchset_alleleseries(observed_guide_reporter_umi_counts_inferred: GeneralMappingInferenceDict, attribute_name: str, contains_guide_surrogate: bool, contains_guide_barcode: bool, contains_guide_umi: bool):
    """Build per-tier observed-allele count Series from a full (retained) mapping result.

    For each whitelist guide, aggregates the observed protospacer/surrogate/barcode
    alleles that mapped to it under the chosen tier, across nine (ambiguity x UMI)
    counting strategies. Required input for downstream mutation profiling.

    Parameters
    ----------
    observed_guide_reporter_umi_counts_inferred
        The per-observation inference dict — i.e. ``result.observed_guide_reporter_umi_counts_inferred``
        from a mapping call made with ``retain_inference_results=True``. A slim
        result (default) passes ``None`` here and raises ``ValueError``.
    attribute_name
        Match tier to extract. Pass a ``MatchTier`` enum member (or its string
        value). Typical: ``MatchTier.PM_SM_BM`` for full-triplet screens.
    contains_guide_surrogate, contains_guide_barcode, contains_guide_umi
        Must match what was configured during mapping (these drive the output
        DataFrame column shape).

    Returns
    -------
    MatchSetWhitelistReporterObservedSequenceCounterSeriesResults
        Dataclass with 9 alleledict + 9 alleleseries_dict + 9 allele_df fields,
        one per (ambiguity, UMI) strategy.

    Raises
    ------
    ValueError
        If called on a slim mapping result (re-run with ``retain_inference_results=True``).
    """
    _require_inference_dict(observed_guide_reporter_umi_counts_inferred, "get_matchset_alleleseries")
    #
    #   DEFINE THE DEFAULTDICTS FOR COUNTING
    #
    ambiguous_ignored_umi_noncollapsed_alleledict : GeneralAlleleDict  = defaultdict(lambda: defaultdict(int))
    ambiguous_ignored_umi_collapsed_alleledict : GeneralAlleleDict  = defaultdict(lambda: defaultdict(int))
    ambiguous_ignored_alleledict : GeneralAlleleDict  = defaultdict(lambda: defaultdict(int))

    ambiguous_accepted_umi_noncollapsed_alleledict : GeneralAlleleDict = defaultdict(lambda: defaultdict(int))
    ambiguous_accepted_umi_collapsed_alleledict : GeneralAlleleDict = defaultdict(lambda: defaultdict(int))
    ambiguous_accepted_alleledict : GeneralAlleleDict = defaultdict(lambda: defaultdict(int))

    ambiguous_spread_umi_noncollapsed_alleledict : GeneralAlleleDict = defaultdict(lambda: defaultdict(int))
    ambiguous_spread_umi_collapsed_alleledict : GeneralAlleleDict = defaultdict(lambda: defaultdict(int))
    ambiguous_spread_alleledict : GeneralAlleleDict = defaultdict(lambda: defaultdict(int))

    #
    # ITERATE THROUGH THE NON-ERROR INFERRED RESULTS AND FILL THE COUNTS
    #
    inferred_value_results: InferenceResult
    for observed_sequence, inferred_value_results in get_non_error_dict(observed_guide_reporter_umi_counts_inferred, attribute_name).items():
        #
        #   Get the relevant attributes
        #
        observed_value_counts: Union[int, CounterType[Optional[str]]] = inferred_value_results.observed_value
        inferred_value_result: CompleteInferenceMatchResult =  inferred_value_results.inferred_value 
            
        match_set_single_inference_match_result : Optional[MatchSetSingleInferenceMatchResult] = getattr(inferred_value_result, attribute_name)
        assert match_set_single_inference_match_result is not None, "match_set_single_inference_match_result should not be none since this is from the non error list. Developer error."

        # §2.4: matches is now a tuple-of-tuples (was Optional[pd.DataFrame]).
        matches: Optional[Tuple[Tuple[Any, ...], ...]] = match_set_single_inference_match_result.value.matches
        if matches:
            num_matches = len(matches)
            for whitelist_sequence_index in matches:
                observed_sequence_index = tuple(observed_sequence)
                if contains_guide_umi:
                    assert isinstance(observed_value_counts, Counter), f"For UMI, expecting observed value is a Counter, but type is {type(observed_value_counts)}"
                    ambiguous_accepted_umi_noncollapsed_alleledict[whitelist_sequence_index][observed_sequence_index] += sum(observed_value_counts.values())
                    ambiguous_accepted_umi_collapsed_alleledict[whitelist_sequence_index][observed_sequence_index] += len(observed_value_counts.values())

                    ambiguous_spread_umi_noncollapsed_alleledict[whitelist_sequence_index][observed_sequence_index] += sum(observed_value_counts.values()) / float(num_matches)
                    ambiguous_spread_umi_collapsed_alleledict[whitelist_sequence_index][observed_sequence_index] += len(observed_value_counts.values()) / float(num_matches)

                    # If there is no ambiguous matches, then add to ambiguous_ignored counter
                    if num_matches == 1:
                        ambiguous_ignored_umi_noncollapsed_alleledict[whitelist_sequence_index][observed_sequence_index] += sum(observed_value_counts.values())
                        ambiguous_ignored_umi_collapsed_alleledict[whitelist_sequence_index][observed_sequence_index] += len(observed_value_counts.values())

                # STANDARD NON-UMI BASED COUNTING
                else:
                    assert isinstance(observed_value_counts, int), f"For non UMI, expecting observed value is an int, but type is {type(observed_value_counts)}"
                    ambiguous_accepted_alleledict[whitelist_sequence_index][observed_sequence_index] += observed_value_counts
                    ambiguous_spread_alleledict[whitelist_sequence_index][observed_sequence_index] += observed_value_counts / float(num_matches)

                    # If there is no ambiguous matches, then add to ambiguous_ignored counter
                    if num_matches == 1:
                        ambiguous_ignored_alleledict[whitelist_sequence_index][observed_sequence_index] += observed_value_counts

                        
    # Helper function that converts defaultdict to series
    def create_dict_counterseries(alleledict: GeneralAlleleDict) -> GeneralAlleleCountSeriesDict: 
        return {whitelist_sequence_key: pd.Series(observed_sequence_counterdict) for whitelist_sequence_key, observed_sequence_counterdict in alleledict.items()}
    
    def create_df_from_dict_counterseries(alleledict_counterseries: GeneralAlleleCountSeriesDict, contains_guide_surrogate: bool, contains_guide_barcode: bool):
        whitelist_reporter_columns = ["protospacer"] # For dynamically creating immutable tuple from reporter sequences based on if provided
                
        if contains_guide_surrogate:
            whitelist_reporter_columns.append("surrogate")
        if contains_guide_barcode:
            whitelist_reporter_columns.append("barcode")


        if alleledict_counterseries is not None:
            observed_reporter_df_counts_list = []
            for whitelist_reporter_tuple, observed_reporter_series_counts in alleledict_counterseries.items():
                observed_reporter_df_counts = pd.DataFrame(observed_reporter_series_counts)
                observed_reporter_df_counts.columns = ["counts"]
                observed_reporter_df_counts_index = observed_reporter_df_counts.index.to_frame()
                observed_reporter_df_counts_index.columns = ["observed_" + label for label in whitelist_reporter_columns]
                observed_reporter_df_counts = pd.concat([observed_reporter_df_counts, observed_reporter_df_counts_index], axis=1)
                observed_reporter_df_counts[["true_" + label for label in whitelist_reporter_columns]] = whitelist_reporter_tuple 
                observed_reporter_df_counts_list.append(observed_reporter_df_counts)
            observed_reporter_df_counts_total = pd.concat(observed_reporter_df_counts_list)
            return observed_reporter_df_counts_total
        else:
            return None
    #
    #   CONVERT THE COUNT DICTS INTO PANDAS SERIES, since this is a more ideal structure.
    #
    match_set_whitelist_reporter_observed_sequence_counter_series_results = MatchSetWhitelistReporterObservedSequenceCounterSeriesResults()
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_umi_noncollapsed_alleleseries_dict = create_dict_counterseries(ambiguous_ignored_umi_noncollapsed_alleledict)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_umi_collapsed_alleleseries_dict = create_dict_counterseries(ambiguous_ignored_umi_collapsed_alleledict)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_alleleseries_dict = create_dict_counterseries(ambiguous_ignored_alleledict)

    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_umi_noncollapsed_alleleseries_dict = create_dict_counterseries(ambiguous_accepted_umi_noncollapsed_alleledict)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_umi_collapsed_alleleseries_dict = create_dict_counterseries(ambiguous_accepted_umi_collapsed_alleledict)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_alleleseries_dict = create_dict_counterseries(ambiguous_accepted_alleledict)

    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_umi_noncollapsed_alleleseries_dict = create_dict_counterseries(ambiguous_spread_umi_noncollapsed_alleledict)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_umi_collapsed_alleleseries_dict = create_dict_counterseries(ambiguous_spread_umi_collapsed_alleledict)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_alleleseries_dict = create_dict_counterseries(ambiguous_spread_alleledict)

    #
    #    Create the dataframes based on the countdicts
    #
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_umi_noncollapsed_allele_df = create_df_from_dict_counterseries(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_umi_noncollapsed_alleleseries_dict, contains_guide_surrogate=contains_guide_surrogate, contains_guide_barcode=contains_guide_barcode)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_umi_collapsed_allele_df = create_df_from_dict_counterseries(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_umi_collapsed_alleleseries_dict, contains_guide_surrogate=contains_guide_surrogate, contains_guide_barcode=contains_guide_barcode)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_allele_df = create_df_from_dict_counterseries(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_alleleseries_dict, contains_guide_surrogate=contains_guide_surrogate, contains_guide_barcode=contains_guide_barcode)

    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_umi_noncollapsed_allele_df = create_df_from_dict_counterseries(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_umi_noncollapsed_alleleseries_dict, contains_guide_surrogate=contains_guide_surrogate, contains_guide_barcode=contains_guide_barcode)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_umi_collapsed_allele_df = create_df_from_dict_counterseries(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_umi_collapsed_alleleseries_dict, contains_guide_surrogate=contains_guide_surrogate, contains_guide_barcode=contains_guide_barcode)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_allele_df = create_df_from_dict_counterseries(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_alleleseries_dict, contains_guide_surrogate=contains_guide_surrogate, contains_guide_barcode=contains_guide_barcode)

    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_umi_noncollapsed_allele_df = create_df_from_dict_counterseries(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_umi_noncollapsed_alleleseries_dict, contains_guide_surrogate=contains_guide_surrogate, contains_guide_barcode=contains_guide_barcode)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_umi_collapsed_allele_df = create_df_from_dict_counterseries(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_umi_collapsed_alleleseries_dict, contains_guide_surrogate=contains_guide_surrogate, contains_guide_barcode=contains_guide_barcode)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_allele_df = create_df_from_dict_counterseries(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_alleleseries_dict, contains_guide_surrogate=contains_guide_surrogate, contains_guide_barcode=contains_guide_barcode)
        
    return match_set_whitelist_reporter_observed_sequence_counter_series_results














from copy import deepcopy
import numpy as np
# TODO: Move to package

def pad_series(series, possible_max = 0):
    max_len = max(series.apply(len).max(), possible_max)
    return series.apply(lambda item: item.ljust(max_len, 'X')), max_len

def pad_sequence(sequence, sequence_max=0):
    return sequence.ljust(sequence_max, "X")

# Copying function from crispr-millipede package - don't want import entire package for a single function
def get_substitution_encoding(true_sequence, observed_sequence, skip_index=0):
    assert len(true_sequence) == len(observed_sequence), f"Sequences not same length. True={true_sequence}, Observed={observed_sequence}" # Ensure the aligned sequence (from allele table) is equal size to the reference sequence
    
    nucleotides = ["A","C","T","G","N","-", "X"] # List of possible nucleotides
    encodings_per_position = []
    mismatch_mappings_per_position = []
    for index in range(0, len(true_sequence)): # Iterate through each base and check for substitution
        # TODO Ensure sequences are uppercase
        nucleotides_mm = nucleotides[:]
        nucleotides_mm.remove(true_sequence[index])
        mm_encoding = pd.Series(np.repeat(0, len(nucleotides)-1))
        if observed_sequence[index] == true_sequence[index]: # If the observed sequence is same as reference
            pass
        else:
            mm_index = nucleotides_mm.index(observed_sequence[index])
            mm_encoding[mm_index] = 1
        mismatch_mappings_per_position.append(nucleotides_mm)
        encodings_per_position.append(mm_encoding)

    encodings_per_position_df = pd.DataFrame(encodings_per_position).T
    mismatch_mappings_per_position_df = pd.DataFrame(mismatch_mappings_per_position).T

    encodings_per_position_df.columns = list(true_sequence)
    mismatch_mappings_per_position_df.columns = list(true_sequence)

    mismatch_mappings_per_position_POS_list = np.arange(mismatch_mappings_per_position_df.shape[1]).repeat(mismatch_mappings_per_position_df.shape[0])
    mismatch_mappings_per_position_REF_list = np.asarray(list(true_sequence)).repeat(mismatch_mappings_per_position_df.shape[0]).astype(np.object_)
    mismatch_mappings_per_position_ALT_list = mismatch_mappings_per_position_df.T.values.flatten()
    mismatch_mappings_per_position_full_list = mismatch_mappings_per_position_POS_list.astype(np.str_).astype(object)+mismatch_mappings_per_position_REF_list + np.repeat(">", len(mismatch_mappings_per_position_REF_list)) + mismatch_mappings_per_position_ALT_list
    encodings_per_position_list = encodings_per_position_df.T.values.flatten()

    
    # Encodings per position DF, mismatch mappings per position DF, encodings per position flattened, mismatch mappings per position flattened, mismatch mapping position in flattened list, mismatch mapping ref in flattened list, mismatch mapping alt in flattened list, all substitutions made
    
    index = pd.MultiIndex.from_tuples(zip(mismatch_mappings_per_position_full_list, mismatch_mappings_per_position_POS_list, mismatch_mappings_per_position_REF_list, mismatch_mappings_per_position_ALT_list), names=["FullChange", "Position","Ref", "Alt"])
    
    assert len(encodings_per_position_list) == len(index), f"{len(encodings_per_position_list)} == {len(index)}"
    encodings_per_position_series = pd.Series(encodings_per_position_list, index = index, name="encoding")
    return encodings_per_position_series


def determine_mutations_in_sequence(true_sequence, observed_sequence):
    mutation_array = []
    lowest_len = min(len(true_sequence), len(observed_sequence))
    for i in range(lowest_len):
        if true_sequence[i] != observed_sequence[i]:
            preceding_nt_context = None
            succeeding_nt_context = None
            ref_nt = true_sequence[i]
            alt_nt = observed_sequence[i]
            if i == 0:
                preceding_nt_context = "_"
                succeeding_nt_context = true_sequence[i+1]
            elif i == len(true_sequence)-1:
                preceding_nt_context = true_sequence[i-1]
                succeeding_nt_context = "_"
            else:
                preceding_nt_context = true_sequence[i-1]
                succeeding_nt_context = true_sequence[i+1]
                
            mutation_series = pd.Series((preceding_nt_context, ref_nt, alt_nt,succeeding_nt_context, i), index=["preceding_nt_context", "ref_nt", "alt_nt", "succeeding_nt_context", "position"])
            mutation_array.append(mutation_series)
    observed_sequence_mutation_df = pd.DataFrame(mutation_array)
    return observed_sequence_mutation_df


from multiprocessing import Pool as _Pool

# ---------------------------------------------------------------------------
# PERF: parallelisation helpers for get_mutation_profile's per-guide loop.
# The per-guide work is embarrassingly parallel across whitelist reporters. The (large)
# alleleseries dict + the two component flags are shared with worker processes via fork
# copy-on-write through these module globals (set by _gmp_set_context just before the Pool is
# created), so only the small per-guide iterrows() row is pickled to each worker. The caller
# accumulates results in input order, so the output is identical to the original serial loop.
# ---------------------------------------------------------------------------
_GMP_ALLELESERIES = None
_GMP_CONTAINS_SURROGATE = False
_GMP_CONTAINS_BARCODE = False
_GMP_COMPUTE_UNLINKED = True   # PERF: skip the (often-discarded) unlinked tables when False
_GMP_MAX_READS = None          # OPT: cap each guide to its top-N reads (most-abundant alleles) if set
_GMP_MAX_ALLELES = None        # OPT: cap each guide to its top-N distinct alleles if set


def _gmp_set_context(alleleseries, contains_guide_surrogate, contains_guide_barcode,
                     compute_unlinked=True, max_reads_per_guide=None, max_alleles_per_guide=None):
    global _GMP_ALLELESERIES, _GMP_CONTAINS_SURROGATE, _GMP_CONTAINS_BARCODE
    global _GMP_COMPUTE_UNLINKED, _GMP_MAX_READS, _GMP_MAX_ALLELES
    _GMP_ALLELESERIES = alleleseries
    _GMP_CONTAINS_SURROGATE = contains_guide_surrogate
    _GMP_CONTAINS_BARCODE = contains_guide_barcode
    _GMP_COMPUTE_UNLINKED = compute_unlinked
    _GMP_MAX_READS = max_reads_per_guide
    _GMP_MAX_ALLELES = max_alleles_per_guide


# Nucleotide order used by get_substitution_encoding (the 6 alts per position = this list minus the ref base).
_GMP_NUC = ["A", "C", "T", "G", "N", "-", "X"]


def get_substitution_encoding_batch(true_sequence, observed_sequences):
    """Vectorized, bit-identical batch equivalent of `get_substitution_encoding` applied to M aligned observed
    sequences at once. Returns an M x 6L DataFrame (RangeIndex rows) whose columns are the MultiIndex
    (FullChange, Position, Ref, Alt) in the SAME order the per-allele function produces (position-major; the 6
    alts per position = `_GMP_NUC` minus the reference base, order preserved); value = 1 iff observed[m][i]==alt
    else 0. The caller prepends the SequenceType level and adds the count column."""
    L = len(true_sequence); M = len(observed_sequences)
    cols = []
    for i in range(L):
        ref = true_sequence[i]
        for alt in [n for n in _GMP_NUC if n != ref]:
            cols.append((f"{i}{ref}>{alt}", i, ref, alt))
    colindex = pd.MultiIndex.from_tuples(cols, names=["FullChange", "Position", "Ref", "Alt"])
    if M == 0:
        return pd.DataFrame(np.empty((0, len(cols)), dtype=np.int64), columns=colindex)
    obs = np.array([list(s) for s in observed_sequences])          # M x L char array
    blocks = []
    for i in range(L):
        ref = true_sequence[i]; col_i = obs[:, i]
        blocks.append(np.stack([(col_i == alt) for alt in [n for n in _GMP_NUC if n != ref]], axis=1).astype(np.int64))
    return pd.DataFrame(np.concatenate(blocks, axis=1), columns=colindex)


def _gmp_process_guide(whitelist_reporter_sequence):
    """Process ONE whitelist reporter row (an iterrows() (index, Series) tuple), returning
    (whitelist_reporter_tuple, linked_mutations_df, protospacer_unlinked_list,
     surrogate_unlinked_list, barcode_unlinked_list) or None if the guide has no observed
    sequences. Verbatim extraction of the original per-guide loop body; reads the shared
    alleleseries + flags from module globals so it is picklable for multiprocessing.Pool."""
    alleleseries = _GMP_ALLELESERIES
    contains_guide_surrogate = _GMP_CONTAINS_SURROGATE
    contains_guide_barcode = _GMP_CONTAINS_BARCODE
    compute_unlinked = _GMP_COMPUTE_UNLINKED
    max_reads = _GMP_MAX_READS
    max_alleles = _GMP_MAX_ALLELES

    protospacer_unlinked_list = []
    surrogate_unlinked_list = []
    barcode_unlinked_list = []

    whitelist_sequence_pretuple_list = []
    whitelist_protospacer_sequence = whitelist_reporter_sequence[1]["protospacer"][:]
    whitelist_sequence_pretuple_list.append(whitelist_protospacer_sequence)
    if contains_guide_surrogate:
        whitelist_surrogate_sequence = whitelist_reporter_sequence[1]["surrogate"][:]
        whitelist_sequence_pretuple_list.append(whitelist_surrogate_sequence)
    if contains_guide_barcode:
        whitelist_barcode_sequence = whitelist_reporter_sequence[1]["barcode"][:]
        whitelist_sequence_pretuple_list.append(whitelist_barcode_sequence)
    whitelist_reporter_tuple = tuple(whitelist_sequence_pretuple_list)

    try:
        observed_sequences_df = alleleseries[whitelist_reporter_tuple]
        observed_sequences_df = observed_sequences_df.copy()  # avoid mutating the shared alleleseries' index
        # OPT: top-N reads / alleles cap (keep the most-abundant alleles; preserve their original ordering)
        if (max_alleles is not None or max_reads is not None) and len(observed_sequences_df) > 0:
            ordered = observed_sequences_df.sort_values(ascending=False, kind="mergesort")
            if max_alleles is not None and len(ordered) > int(max_alleles):
                ordered = ordered.iloc[:int(max_alleles)]
            if max_reads is not None and float(ordered.sum()) > float(max_reads):
                csum = np.cumsum(ordered.to_numpy())
                keep = int(np.searchsorted(csum, float(max_reads), side="left")) + 1
                ordered = ordered.iloc[:keep]
            if len(ordered) < len(observed_sequences_df):
                observed_sequences_df = observed_sequences_df[observed_sequences_df.index.isin(set(ordered.index))]
        whitelist_reporter_sequence_copy = deepcopy(whitelist_reporter_sequence)
        padded_observed_sequences_results = []
        for level in range(observed_sequences_df.index.nlevels):
            series_padding_results = pad_series(observed_sequences_df.index.get_level_values(level).to_series(), possible_max=len(whitelist_reporter_sequence_copy[1].iloc[level]))
            padded_observed_sequences_results.append(series_padding_results[0].to_list())
            whitelist_reporter_sequence_copy[1].iloc[level] = pad_sequence(whitelist_reporter_sequence_copy[1].iloc[level], sequence_max=series_padding_results[1])
        padded_observed_sequence_results_tuples = [tuple(padded_observed_sequence_list) for padded_observed_sequence_list in zip(*padded_observed_sequences_results)]
        observed_sequences_df.index = pd.MultiIndex.from_tuples(padded_observed_sequence_results_tuples)

        # ---- LINKED (vectorized, bit-identical): per component batch-encode, prepend SequenceType, add count ----
        obs_index = list(observed_sequences_df.index)
        counts = observed_sequences_df.to_numpy()
        comp_defs = [("protospacer", 0)]
        if contains_guide_surrogate: comp_defs.append(("surrogate", 1))
        if contains_guide_barcode: comp_defs.append(("barcode", 2))
        parts = []
        for name, pos in comp_defs:
            enc = get_substitution_encoding_batch(whitelist_reporter_sequence_copy[1][name], [o[pos] for o in obs_index])
            enc.columns = pd.MultiIndex.from_tuples([(name,) + tuple(c) for c in enc.columns],
                                                    names=["SequenceType", "FullChange", "Position", "Ref", "Alt"])
            parts.append(enc)
        count_df = pd.DataFrame(counts.reshape(-1, 1),
                                columns=pd.MultiIndex.from_tuples([("count", "", "", "", "")],
                                                                 names=["SequenceType", "FullChange", "Position", "Ref", "Alt"]))
        linked_mutations_df = pd.concat(parts + [count_df], axis=1)
        linked_mutations_df.index = observed_sequences_df.index

        # ---- UNLINKED (only when requested; some callers, e.g. the linked-only sidecar path, discard these) ----
        if compute_unlinked:
            for observed_sequences, count in observed_sequences_df.items():
                observed_protospacer_unlinked_mutations_df = determine_mutations_in_sequence(true_sequence=whitelist_reporter_sequence_copy[1]["protospacer"], observed_sequence=observed_sequences[0])
                observed_protospacer_unlinked_mutations_df["count"] = count
                observed_protospacer_unlinked_mutations_df.loc[:, whitelist_reporter_sequence_copy[1].index.values] = whitelist_reporter_sequence_copy[1].values
                protospacer_unlinked_list.append(observed_protospacer_unlinked_mutations_df)
                if contains_guide_surrogate:
                    observed_surrogate_unlinked_mutations_df = determine_mutations_in_sequence(true_sequence=whitelist_reporter_sequence_copy[1]["surrogate"], observed_sequence=observed_sequences[1])
                    observed_surrogate_unlinked_mutations_df["count"] = count
                    observed_surrogate_unlinked_mutations_df.loc[:, whitelist_reporter_sequence_copy[1].index.values] = whitelist_reporter_sequence_copy[1].values
                    surrogate_unlinked_list.append(observed_surrogate_unlinked_mutations_df)
                if contains_guide_barcode:
                    observed_barcode_unlinked_mutations_df = determine_mutations_in_sequence(true_sequence=whitelist_reporter_sequence_copy[1]["barcode"], observed_sequence=observed_sequences[2])
                    observed_barcode_unlinked_mutations_df["count"] = count
                    observed_barcode_unlinked_mutations_df.loc[:, whitelist_reporter_sequence_copy[1].index.values] = whitelist_reporter_sequence_copy[1].values
                    barcode_unlinked_list.append(observed_barcode_unlinked_mutations_df)
        return (whitelist_reporter_tuple, linked_mutations_df, protospacer_unlinked_list, surrogate_unlinked_list, barcode_unlinked_list)
    except KeyError:
        print(f"No observed sequences found for: {whitelist_reporter_tuple}")
        return None


def get_mutation_profile(match_set_whitelist_reporter_observed_sequence_counter_series_results: MatchSetWhitelistReporterObservedSequenceCounterSeriesResults, whitelist_reporter_df: pd.DataFrame, contains_guide_surrogate: bool, contains_guide_barcode: bool, strategies: Optional[List[str]] = None, cores: int = 1, compute_unlinked: bool = True, max_reads_per_guide: Optional[int] = None, max_alleles_per_guide: Optional[int] = None) -> MatchSetWhitelistReporterObservedSequenceMutationProfiles:
    """Compute per-position mutation profiles from allele count series.

    Given the allele Series built by ``get_matchset_alleleseries``, this walks
    each (whitelist, observed_allele) pair and records per-base mutations
    against the whitelist reference, producing both linked (allele-level) and
    unlinked (position-level) mutation tables for all nine counting strategies.

    Parameters
    ----------
    match_set_whitelist_reporter_observed_sequence_counter_series_results
        Return value of ``get_matchset_alleleseries``.
    whitelist_reporter_df
        The same DataFrame that was passed into the mapping call. Used as the
        reference sequence for computing mutations.
    contains_guide_surrogate, contains_guide_barcode
        Must match the mapping configuration.
    strategies
        Optional list of strategy base-names (e.g. ``["ambiguous_accepted_umi_noncollapsed"]``)
        to compute; ``None`` (default) computes all nine, reproducing the original behaviour.
    cores
        Number of worker processes for the per-guide loop within each strategy. ``1`` (default)
        runs serially (identical to the original); ``>1`` parallelises via multiprocessing.

    Returns
    -------
    MatchSetWhitelistReporterObservedSequenceMutationProfiles
        Mutation tables per strategy. Consume via
        ``tally_linked_mutation_count_per_sequence`` for aggregate counters, or
        drive ``visualization.plot_mutation_count_histogram`` /
        ``plot_trinucleotide_mutational_signature`` directly.

    See Also
    --------
    tally_linked_mutation_count_per_sequence
    """
    # Function to generate unlinked mutations for particular count type
    def generate_mutations_results(alleleseries: Optional[GeneralAlleleCountSeriesDict], whitelist_reporter_df: pd.DataFrame, contains_guide_surrogate: bool, contains_guide_barcode: bool) -> Optional[ObservedSequenceMutationProfile]:
        if alleleseries is not None:
            
            linked_mutations_whitelist_reporter_dict = {}
            all_observed_protospacer_unlinked_mutations_df_list = []
            if contains_guide_surrogate:
                all_observed_surrogate_unlinked_mutations_df_list = []
            if contains_guide_barcode:
                all_observed_barcode_unlinked_mutations_df_list = []

            _gmp_set_context(alleleseries, contains_guide_surrogate, contains_guide_barcode,
                             compute_unlinked=compute_unlinked, max_reads_per_guide=max_reads_per_guide,
                             max_alleles_per_guide=max_alleles_per_guide)
            _rows = list(whitelist_reporter_df.iterrows())
            if cores and cores > 1 and len(_rows) > 1:                     # PERF: parallel per-guide loop
                _nproc = min(int(cores), len(_rows))
                with _Pool(_nproc) as _pool:
                    _per_guide_results = _pool.map(_gmp_process_guide, _rows, chunksize=max(1, len(_rows) // (_nproc * 4)))
            else:
                _per_guide_results = [_gmp_process_guide(_row) for _row in _rows]

            for _res in _per_guide_results:                                # accumulate in input order -> identical to serial
                if _res is None:
                    continue
                _wl_tuple, _linked_df, _proto_list, _surr_list, _bc_list = _res
                linked_mutations_whitelist_reporter_dict.update({_wl_tuple: _linked_df})
                all_observed_protospacer_unlinked_mutations_df_list.extend(_proto_list)
                if contains_guide_surrogate:
                    all_observed_surrogate_unlinked_mutations_df_list.extend(_surr_list)
                if contains_guide_barcode:
                    all_observed_barcode_unlinked_mutations_df_list.extend(_bc_list)
                
                
                
            # unlinked tables only when requested (compute_unlinked=False leaves them None); guard empty concat
            all_observed_protospacer_unlinked_mutations_df = (pd.concat(all_observed_protospacer_unlinked_mutations_df_list)
                                                              if (compute_unlinked and all_observed_protospacer_unlinked_mutations_df_list) else None) # Concat list of all whitelist mutations to sequence dataframe
            observed_sequence_mutations = ObservedSequenceMutationProfile(all_observed_protospacer_unlinked_mutations_df=all_observed_protospacer_unlinked_mutations_df, linked_mutations_whitelist_reporter_dict=linked_mutations_whitelist_reporter_dict) # Initialize results object with groupby results

            # Add surrogate result to result object (if provided)
            if contains_guide_surrogate and compute_unlinked and all_observed_surrogate_unlinked_mutations_df_list:
                all_observed_surrogate_unlinked_mutations_df = pd.concat(all_observed_surrogate_unlinked_mutations_df_list)
                observed_sequence_mutations.all_observed_surrogate_unlinked_mutations_df=all_observed_surrogate_unlinked_mutations_df

            # Add barcode result to result object (if provided)
            if contains_guide_barcode and compute_unlinked and all_observed_barcode_unlinked_mutations_df_list:
                all_observed_barcode_unlinked_mutations_df = pd.concat(all_observed_barcode_unlinked_mutations_df_list)
                observed_sequence_mutations.all_observed_barcode_unlinked_mutations_df=all_observed_barcode_unlinked_mutations_df
            
            return observed_sequence_mutations
        else:
            return None
    
    # Generate mutation counts for each count type
    mutations_results = MatchSetWhitelistReporterObservedSequenceMutationProfiles()
    def _want(_name):                                                  # PERF: compute only requested strategies
        return strategies is None or _name in strategies
    if _want("ambiguous_ignored_umi_noncollapsed"):
        print("Generating ambiguous_ignored_umi_noncollapsed_mutations")
        mutations_results.ambiguous_ignored_umi_noncollapsed_mutations = generate_mutations_results(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_umi_noncollapsed_alleleseries_dict, whitelist_reporter_df, contains_guide_surrogate, contains_guide_barcode)
    if _want("ambiguous_ignored_umi_collapsed"):
        print("Generating ambiguous_ignored_umi_collapsed_mutations")
        mutations_results.ambiguous_ignored_umi_collapsed_mutations = generate_mutations_results(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_umi_collapsed_alleleseries_dict, whitelist_reporter_df, contains_guide_surrogate, contains_guide_barcode)
    if _want("ambiguous_ignored"):
        print("Generating ambiguous_ignored_unlinked_mutations")
        mutations_results.ambiguous_ignored_mutations = generate_mutations_results(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_alleleseries_dict, whitelist_reporter_df, contains_guide_surrogate, contains_guide_barcode)

    if _want("ambiguous_accepted_umi_noncollapsed"):
        print("Generating ambiguous_accepted_umi_noncollapsed_mutations")
        mutations_results.ambiguous_accepted_umi_noncollapsed_mutations = generate_mutations_results(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_umi_noncollapsed_alleleseries_dict, whitelist_reporter_df, contains_guide_surrogate, contains_guide_barcode)
    if _want("ambiguous_accepted_umi_collapsed"):
        print("Generating ambiguous_accepted_umi_collapsed_mutations")
        mutations_results.ambiguous_accepted_umi_collapsed_mutations = generate_mutations_results(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_umi_collapsed_alleleseries_dict, whitelist_reporter_df, contains_guide_surrogate, contains_guide_barcode)
    if _want("ambiguous_accepted"):
        print("Generating ambiguous_accepted_mutations")
        mutations_results.ambiguous_accepted_mutations = generate_mutations_results(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_alleleseries_dict, whitelist_reporter_df, contains_guide_surrogate, contains_guide_barcode)

    if _want("ambiguous_spread_umi_noncollapsed"):
        print("Generating ambiguous_spread_umi_noncollapsed_mutations")
        mutations_results.ambiguous_spread_umi_noncollapsed_mutations = generate_mutations_results(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_umi_noncollapsed_alleleseries_dict, whitelist_reporter_df, contains_guide_surrogate, contains_guide_barcode)
    if _want("ambiguous_spread_umi_collapsed"):
        print("Generating ambiguous_spread_umi_collapsed_mutations")
        mutations_results.ambiguous_spread_umi_collapsed_mutations = generate_mutations_results(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_umi_collapsed_alleleseries_dict, whitelist_reporter_df, contains_guide_surrogate, contains_guide_barcode)
    if _want("ambiguous_spread"):
        print("Generating ambiguous_spread_mutations")
        mutations_results.ambiguous_spread_mutations = generate_mutations_results(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_alleleseries_dict, whitelist_reporter_df, contains_guide_surrogate, contains_guide_barcode)
    
    return mutations_results

 
def tally_linked_mutation_count_per_sequence(mutations_results: MatchSetWhitelistReporterObservedSequenceMutationProfiles, contains_guide_surrogate: bool, contains_guide_barcode: bool, count_attribute_name: str="ambiguous_accepted_umi_noncollapsed_mutations") -> LinkedMutationCounters:
    protospacer_total_mutation_counter: CounterType = Counter()

    surrogate_total_mutation_counter: Optional[CounterType] = None
    barcode_total_mutation_counter: Optional[CounterType] = None
    if contains_guide_surrogate:
        surrogate_total_mutation_counter = Counter()
    if contains_guide_barcode:
        barcode_total_mutation_counter = Counter()

    # TODO: THIS IS HARDCODED TO ambiguous_accepted_umi_noncollapsed
        
    mutation_results_typed = getattr(mutations_results, count_attribute_name)
    if mutation_results_typed is not None:
        for linked_mutations_whitelist_reporter_df in mutation_results_typed.linked_mutations_whitelist_reporter_dict.values():
            accetable_SNP_columns = (linked_mutations_whitelist_reporter_df.columns.get_level_values("SequenceType") != "count") & (linked_mutations_whitelist_reporter_df.columns.get_level_values("Ref") != "X") & (linked_mutations_whitelist_reporter_df.columns.get_level_values("Alt") != "X")


            protospacer_columns = (linked_mutations_whitelist_reporter_df.columns.get_level_values("SequenceType") == "protospacer")
            protospacer_mutation_count_per_allele_series = linked_mutations_whitelist_reporter_df.loc[:, protospacer_columns & accetable_SNP_columns].sum(axis=1)
            for allele_index, mutation_count in enumerate(protospacer_mutation_count_per_allele_series):
                protospacer_total_mutation_counter[mutation_count] += linked_mutations_whitelist_reporter_df["count"][allele_index]

            if contains_guide_surrogate:
                surrogate_columns = (linked_mutations_whitelist_reporter_df.columns.get_level_values("SequenceType") == "surrogate")
                surrogate_mutation_count_per_allele_series = linked_mutations_whitelist_reporter_df.loc[:, surrogate_columns & accetable_SNP_columns].sum(axis=1)
                for allele_index, mutation_count in enumerate(surrogate_mutation_count_per_allele_series):
                    surrogate_total_mutation_counter[mutation_count] += linked_mutations_whitelist_reporter_df["count"][allele_index]

            if contains_guide_barcode:
                barcode_columns = (linked_mutations_whitelist_reporter_df.columns.get_level_values("SequenceType") == "barcode")
                barcode_mutation_count_per_allele_series = linked_mutations_whitelist_reporter_df.loc[:, barcode_columns & accetable_SNP_columns].sum(axis=1)
                for allele_index, mutation_count in enumerate(barcode_mutation_count_per_allele_series):
                    barcode_total_mutation_counter[mutation_count] += linked_mutations_whitelist_reporter_df["count"][allele_index]
    
    return LinkedMutationCounters(protospacer_total_mutation_counter=protospacer_total_mutation_counter,
                                surrogate_total_mutation_counter=surrogate_total_mutation_counter,
                                barcode_total_mutation_counter=barcode_total_mutation_counter)
