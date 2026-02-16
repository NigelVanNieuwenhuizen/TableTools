#!/usr/bin/env python3

# TableTools
# Author: Nigel Van Nieuwenhuizen, MSc
# Created: February 2021
# Last Updated: February 2026

# Global imports
import math, random, itertools, operator, statistics, datetime, copy, os, io

# Toolbox classes
class _ListOps(): 
    """A set of general functions for manipulating lists and nested lists."""
    # ============================================================
    # PRIVATE METHODS AND CLASS ATTRIBUTES AND FUNCTION DISPLAY
    # ============================================================
    def __init__(self):
        self._total_funcs = len([f for f in dir(__class__) if not f.startswith('__')])
    
    def __repr__(self):
        """If the toolbox is printed, display a message."""
        return ".listops: a set of general functions for manipulating lists and nested lists."

    def display_toolbox_functions(self):
        """Display a list of all available functions within this toolbox."""
        print(f"Number of {__class__.__name__[1:]} functions: {len([f for f in dir(__class__) if not f.startswith('__')])}")
        for f in [f for f in dir(__class__) if not f.startswith("__")]:
            print(f)

    # ============================================================
    # DATA TYPE AND TRANSFORMATION
    # ============================================================

    def convert_data_type(self,vals,dtype = "integer"):
        """Convert each value in an input list ('vals') to the specified data type ('dtype') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        dtype (string): the data type to convert values to. May be specified as "integer" or "float" or "string" or "Boolean" for integer, float, string, or Boolean types, respectively."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if dtype not in {"integer", "float", "string", "boolean"}:
            raise ValueError("'dtype' must be one of: 'integer', 'float', 'string', 'Boolean'.")

        if dtype == "integer":
            return [int(v) for v in vals]
        elif dtype == "float":
            return [float(v) for v in vals]
        elif dtype == "string":
            return [str(v) for v in vals]
        elif dtype == "Boolean":
            return [bool(v) for v in vals]
      
    def is_type(self, vals, target_type="int"):
        """Check the type of each element in a list ('vals') against a target type and return a list of Boolean values.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        target_type (string): the type to check against. May be specified as 'int', 'float', 'string', or 'bool'. """

        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if target_type not in {"bool", "int", "float", "string"}:
            raise ValueError(f"Unknown target_type '{target_type}'. Must be one of 'bool', 'int', 'float', 'string'.")

        type_map = {
            "bool": bool,
            "int": int,
            "float": float,
            "string": str
        }

        check_type = type_map[target_type]
        return [isinstance(v, check_type) for v in vals]
    
    def categorical_encoding(self,vals,method):
        """Encode a list of categorical values ('vals') with the one-hot, dummy, or label encoding method and return a new list.
            
        Parameters:
            
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
            
        method (string): the method of encoding. Possible methods include 'one-hot' or 'dummy' or 'label' encoding."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if method not in {"one-hot", "dummy", "label"}:
            raise ValueError("'method' must be one of: 'one-hot', 'dummy', 'label'.")

        output = []
        lv = len(vals)
        if method == "one-hot":
            for v in range(lv):
                encoded = ["0"] * lv
                encoded[v] = "1"
                output.append("".join(encoded))
        elif method == "dummy":
            for v in range(lv):
                encoded = ["0"] * (lv - 1)
                if v < lv - 1:
                    encoded[v] = "1"
                output.append("".join(encoded))
        elif method == "label":
            categories = {val: idx for idx, val in enumerate(sorted(set(vals)))}
            output = [categories[val] for val in vals]
        return output

    # ============================================================
    # DUPLICATION AND REPLICATION
    # ============================================================

    def duplicate_value(self,vals,dup_value,duplicates):
        """Duplicate all instances of a specified value ('dup_val') from the input list ('vals') a specified number of times ('duplicates') and return a new list. 
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        dup_value (integer, float, string): the value to be duplicated. 
        
        duplicates (integer): the number of times the specified value will appear in the output list, for each instance of the value in the input list."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not isinstance(duplicates, int) or duplicates < 1:
            raise ValueError("'duplicates' must be an integer greater than or equal to 1.")
        output = []
        append = output.append
        for v in vals:
            if v == dup_value:
                for _ in range(duplicates):
                    append(v)
            else:
                append(v)
        return output

    def duplicate_all(self,vals,duplicates):
        """Duplicate all values in the input list ('vals') a specified number of times ('duplicates') and return a new list.
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        duplicates (integer): the number of times each value in the input list will appear in the output list."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not isinstance(duplicates, int) or duplicates < 1:
            raise ValueError("'duplicates' must be an integer greater than or equal to 1.")
        return [v for v in vals for _ in range(duplicates)]

    def remove_duplicates_unordered(self,vals):
        """Remove duplicate values from the input list ('vals') and return a new list. Order of list elements is not preserved.
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        return list(set(vals))

    def remove_duplicates_ordered(self,vals):
        """Remove duplicate values from the input list ('vals') and return a new list. Order of list elements is preserved.
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        seen = set()
        output = []
        for v in vals:
            if v not in seen:
                seen.add(v)
                output.append(v)
        return output
    
    def find_duplicates(self,vals):
        """Find the duplicate values of an input list ('vals') and return a new list.
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        seen = set()
        dups = set()
        for v in vals:
            if v in seen:
                dups.add(v)
            else:
                seen.add(v)
        return list(dups)
    
    def has_duplicates(self,vals):
        """Determine if an input list ('vals') has duplicate values and return True or False.
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        return len(vals) != len(set(vals))

    # ============================================================
    # REMOVAL AND REPLACEMENT
    # ============================================================

    def remove_value(self,vals,to_remove):
        """Remove all instances of a specified value ('to_remove') from the input list ('vals') and return a new list.
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        to_remove (integer, float, string): the value to be removed."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        return [v for v in vals if v != to_remove]

    def remove_nth(self,vals,nth):
        """Remove every nth value ('nth') from the input list ('vals') and return a new list.
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        nth (integer): the nth value to remove. Must be greater than 1."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not isinstance(nth, int) or nth <= 1:
            raise ValueError("'nth' must be an integer greater than 1.")
        v = vals[:]  # copy to avoid mutating input
        del v[nth-1::nth]
        return v
    
    def keep_nth(self, vals, nth):
        """Keep every nth value ('nth') from the input list ('vals') and return a new list.
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered
        by the operation.
        
        nth (integer): the nth value to keep. Must be greater than 1."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not isinstance(nth, int) or nth <= 1:
            raise ValueError("'nth' must be an integer greater than 1.")
        return vals[nth-1::nth]

    def replace_value(self,vals,old_val,new_val):
        """Replace all instances of a specified value ('old_val') with a new value ('new_val') from the input list ('vals') and return a new list.
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        old_val (integer, float, string): the old value to be replaced.
        
        new_val (integer, float, string): the new value that will replace the old value."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        return [v if v != old_val else new_val for v in vals]

    # ============================================================
    # STRUCTURAL OPERATIONS
    # ============================================================    

    def shuffle_list(self,vals):
        """Randomize the order of the input list ('vals') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        return random.sample(vals, len(vals))
  
    def shorten_list(self,vals,new_length,position="end"):
        """Reduce the length of the input list ('vals') to a specified length ('new_length') starting from the specified position ('position') and return a new list.
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        new_length (integer): the new length of the output list. Values in the input list beyond the specified length will be exluded from the new list.
        
        position (string): flag to indicate if values should be removed from the start or the end of the input list. Possible values are 'start' and 'end'"""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not isinstance(new_length, int) or new_length < 0:
            raise ValueError("'new_length' must be a non-negative integer.")
        if new_length > len(vals):
            raise ValueError("'new_length' cannot be greater than the length of 'vals'.")
        if position not in {"start", "end"}:
            raise ValueError("'position' must be either 'start' or 'end'.")
        if position == "end":
            return vals[:new_length]
        else:  # position == "start"
            dex = len(vals) - new_length
            return vals[dex:]
    
    def pad_list(self,vals,new_len,pad_val=""):
        """Add a specified value ('pad_val') to the end of the input list ('vals') and return a new list.
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        new_len (integer): the new length of the output list. Must be longer than the input list.
        
        pad_val (integer, float, string): the value to add to the end of the input list until the length of the input list is equal to the second list."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not isinstance(new_len, int) or new_len <= len(vals):
            raise ValueError("'new_len' must be an integer greater than the length of 'vals'.")
        output = vals[:]  # copy to avoid mutating input
        len_vals = len(vals)
        return output + [pad_val for _ in range(new_len - len_vals)]

    def extract_evenly_spaced_values(self,vals, k):
        """Extract k evenly-spaced values from an input list of values ('vals') and return a list of values. The first and last value of the input list will be included
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        k (integer): the number of values to extract from the input list, including the first and last."""
        if not isinstance(vals, list) or not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not isinstance(k, int):
            raise ValueError("'k' must be an integer.")
        if k < 2:
            raise ValueError("'k' must be at least 2.")
        if k > len(vals):
            return vals[:]
        
        indices = [int(round(i * (len(vals) - 1) / (k - 1)))for i in range(k)]
        return [vals[i] for i in indices]

    def extract_evenly_spaced_indices(self,vals, k):
        """Extract k evenly-spaced value indices from an input list of values ('vals') and return a list of indices. The first and last index of the input list will be included
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        k (integer): the number of indices to extract from the input list, including the first and last."""
        if not isinstance(vals, list) or not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not isinstance(k, int):
            raise ValueError("'k' must be an integer.")
        if k < 2:
            raise ValueError("'k' must be at least 2.")
        if k >= len(vals):
            return list(range(len(vals)))
        
        return [int(round(i * (len(vals) - 1) / (k - 1)))for i in range(k)]


    # ============================================================
    # GROUPING OPERATIONS
    # ============================================================
   
    def group_consecutive_duplicates_ordered(self,vals):
        """Group consecutive duplicates in an input list ('vals') and return a new nested list. The order the duplicates appear in the input list is preserved.
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        return [list(j) for _, j in itertools.groupby(vals)]
    
    def group_consecutive_duplicates_unordered(self,vals):
        """Group consecutive duplicates in an input list ('vals') and return a new nested list. The order the duplicates appear in the input list is not preserved.
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        from collections import Counter
        counts = Counter(vals)
        return [[k] * v for k, v in counts.items()]
  
    def group_by_equal_count(self,vals,sub_lists):
        """Convert an input one-dimensional list ('vals') into a nested list with a specified number of sub-lists ('sub_lists') and return a new list. Equivalent to an equal count/quantile grouping (i.e. each sub-list will have an equal number of values).

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        sub_lists (integer): the number of sub-lists to create from the input list. If values of the input list do not divide evenly between the number of sub-lists specified, the final sub-list in the nested list will contain fewer values than the others."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not isinstance(sub_lists, int) or sub_lists < 1:
            raise ValueError("'sub_lists' must be an integer greater than or equal to 1.")
        sub_count = math.ceil(len(vals) / sub_lists)
        output = [vals[i:i+sub_count] for i in range(0, len(vals), sub_count)]
        return output

    def group_by_equal_interval(self,vals,intervals):
        """Convert an input one-dimensional list ('vals') into a nested list with a specified number of sub-lists ('intervals') and return a new list. Equivalent to an equal interval grouping (i.e. sub-lists will have an equal range of values, but may be of different lengths).

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        intervals (integer): the number of approximately equal-range intervals to create from the input list."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not isinstance(intervals, int) or intervals < 1:
            raise ValueError("'intervals' must be an integer greater than or equal to 1.")
        maxv = max(vals)
        minv = min(vals)
        interval_width = (maxv - minv + 1e-9) / intervals  # tiny epsilon to avoid edge-case overflow
        interval_data = [[] for _ in range(intervals)]
        for v in vals:
            v_int = math.floor((v - minv) / interval_width)
            if v_int == intervals:  # edge case: max value falls exactly on upper bound
                v_int -= 1
            interval_data[v_int].append(v)
        return interval_data

    def group_by_sublist_len(self,vals,nested):
        """Group values in an input one-dimensional list ('vals') into a nested list whose sublists match the lengths of sublists in an input nested list ('nested') and return a new nested list.
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        nested (list): the nested list to match the output sublist length to. Total number of values in all sublists must be equal to the length of the input list, which can be checked with the total_vals_in_nested function."""
        ln = self.len_sublists(nested)  # lengths of each sublist
        output = []
        pos = 0
        for length in ln:
            output.append(vals[pos:pos+length])
            pos += length
        return output

    # ============================================================
    # NESTED LIST OPERATIONS
    # ============================================================

    def flatten_nested_list(self,vals):
        """Convert an input nested list ('vals') into a one-dimensional list and return a new list.

        Parameters:

        vals (list): the list of lists the operation will be performed on. The input list will not be altered by the operation."""
        return list(itertools.chain.from_iterable(vals))
    
    def nest_list(self, vals, sub_len):
        """Convert an input one-dimensional list ('vals') into a nested list of sub-lists, each of a specified length ('sub_len'), and return a new list.
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        sub_len (integer): the length of each sub-list. Must be greater than 0. The final sub-list may contain fewer values if the length of 'vals' is not evenly divisible by 'sub_len'."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not isinstance(sub_len, int) or sub_len <= 0:
            raise ValueError("'sub_len' must be an integer greater than 0.")
        return [vals[i:i+sub_len] for i in range(0, len(vals), sub_len)]

    def total_vals_in_nested(self,vals):
        """Calculate and return the total number of values in all sublists of an input nested list ('vals')."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, list) for v in vals):
            raise ValueError("All elements of 'vals' must be lists.")
        return sum(len(v) for v in vals)
 
    def len_sublists(self,vals):
        """Calculate the length of each sub-list in the input nested list ('vals') and return a new list of lengths.
        
        Parameters:
        
        vals (list):the list of values the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, list) for v in vals):
            raise ValueError("All elements of 'vals' must be lists.")
        return [len(v) for v in vals]
    
    def max_nested_depth(self, vals):
        """Calculate the maximum depth of nesting in an input list ('vals') and return the depth as an integer.
        
        Parameters:
        
        vals (list): the list to be analyzed. The input list will not be altered by the operation."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        def depth(v):
            if isinstance(v, list) and v:
                return 1 + max(depth(i) for i in v)
            elif isinstance(v, list):
                return 1  # empty sub-list still counts as depth 1
            else:
                return 0

        return depth(vals)

    # ============================================================
    # SET AND BOOLEAN OPERATIONS
    # ============================================================ 
        
    def bool_op(self, vals, vals2=None, op="and", numeric=True):
        """Apply a Boolean operation element-wise to one or two lists of Boolean values.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        vals2 (list): the list of values to compare against the list of operation. The input list will not be altered by the operation.

        op (string): the Boolean operation to apply. May be specified as 'and', 'or', 'xor', 'not', 'eq', or 'ne'.
        numeric (bool): flag to indicate if the output values will be represented numerically (1 and 0) or as Boolean (True and False). Default is True."""

        # Validation
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if op in {"and", "or", "xor", "eq", "ne"}:
            if vals2 is None:
                raise ValueError(f"'vals2' must be provided for '{op}' operation.")
            if not isinstance(vals2, list):
                raise ValueError("'vals2' must be a list.")
            if len(vals) != len(vals2):
                raise ValueError("Input lists must be the same length.")

        # Normalize inputs to bool
        vals_bool = [bool(v) for v in vals]
        vals2_bool = [bool(v) for v in vals2] if vals2 is not None else None

        # Apply operation
        if op == "and":
            result = [a and b for a, b in zip(vals_bool, vals2_bool)]
        elif op == "or":
            result = [a or b for a, b in zip(vals_bool, vals2_bool)]
        elif op == "xor":
            result = [a ^ b for a, b in zip(vals_bool, vals2_bool)]
        elif op == "eq":
            result = [a == b for a, b in zip(vals_bool, vals2_bool)]
        elif op == "ne":
            result = [a != b for a, b in zip(vals_bool, vals2_bool)]
        elif op == "not":
            result = [not a for a in vals_bool]
        else:
            raise ValueError(f"Unknown operation '{op}'.")

        # Return numeric or Boolean
        return [int(r) if numeric else r for r in result]
    
    def set_operation(self, vals1, vals2, method="union"):
        """Perform a set operation on two input lists ('vals1' and 'vals2') and return a new list.
        
        Parameters:
        
        vals1 (list): the first list of values the operation will be performed on. The input list will not be altered by the operation.
        
        vals2 (list): the second list of values the operation will be performed on. The input list will not be altered by the operation.
        
        method (string): the set operation to perform. Possible values are "union", "intersection", "difference", and "sym_diff"."""
        if not isinstance(vals1, list) or not isinstance(vals2, list):
            raise ValueError("'vals1' and 'vals2' must both be lists.")
        if method not in {"union", "intersection", "difference", "sym_diff"}:
            raise ValueError("'method' must be one of: 'union', 'intersection', 'difference', 'sym_diff'.")
        set1, set2 = set(vals1), set(vals2)
        if method == "union":
            return list(set1 | set2)
        elif method == "intersection":
            return list(set1 & set2)
        elif method == "difference":
            return list(set1 - set2)
        elif method == "sym_diff":
            return list(set1 ^ set2)

class _Generate():
    """A set of functions for generating new lists of different data types and distributions."""
    # ============================================================
    # PRIVATE METHODS AND CLASS ATTRIBUTES AND FUNCTION DISPLAY
    # ============================================================
    def __init__(self):
        self._total_funcs = len([f for f in dir(__class__) if not f.startswith('__')])
    
    def __repr__(self):
        """If the toolbox is printed, display a message."""
        return ".generate: a set of functions for generating new lists of different data types and distributions."

    def display_toolbox_functions(self):
        """Display a list of all available functions within this toolbox."""
        print(f"Number of {__class__.__name__[1:]} functions: {len([f for f in dir(__class__) if not f.startswith('__')])}")
        for f in [f for f in dir(__class__) if not f.startswith("__")]:
            print(f)

    # ============================================================
    # DETERMINISTIC SEQUENCES
    # ============================================================

    def constant(self,val,length):
        """Return a list of a single specified value ('val') of the specified length ('length').
        
        Parameters:
        
        val (integer, float, string, list, object): the constant value or object of the output list.
        
        length (integer): the number of values or objects in the output list."""
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        return [val for _ in range(length)]

    def arithmetic_seq(self,start,step,length):
        """Return a list of values of the specified length ('length') in which the specified difference between subsequent values ('step') is constant (by addition or subtraction).

        Parameters:

        start (integer, float): the first value of the output data range.

        step (integer, float): the difference between subsequent values. May be positive (for addition) or negative (for subtraction).

        length (integer): the number of values in the output list."""
        if not isinstance(start, (int, float)):
            raise ValueError("'start' must be an integer or float.")
        if not isinstance(step, (int, float)):
            raise ValueError("'step' must be an integer or float.")
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        output = [start]
        for i in range(1, length):
            output.append(output[i - 1] + step)
        return output
    
    def varied_arithmetic_seq(self,start,step,max_val,length):
        """Return a list of values of the specified length ('length') in which the specified difference between subsequent values ('step') is allowed to vary within a specified limit ('max_val').

        Parameters:

        start (integer, float): the first value of the output data range.

        step (integer, float): the difference between subsequent values. May be positive (for addition) or negative (for subtraction).

        max_val (integer): the maximum value by which the step value is allowed to increase or decrease. Smaller values result in less variation.

        length (integer): the number of values in the output list."""
        if not isinstance(start, (int, float)):
            raise ValueError("'start' must be an integer or float.")
        if not isinstance(step, (int, float)):
            raise ValueError("'step' must be an integer or float.")
        if not isinstance(max_val, int) or max_val < 0:
            raise ValueError("'max_val' must be a non-negative integer.")
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        output = [start]
        for i in range(1, length):
            varied_step = random.randint(int(step - max_val), int(step + max_val))
            output.append(output[i - 1] + varied_step)
        return output
    
    def harmonic_seq(self,start,step,length):
        """Return a list of values of the specified length ('length') in which the specified difference between subsequent values ('step') is the reciprocal of an arithmetic sequence.

        Parameters:

        start (integer, float): the first value of the output data range.

        step (integer, float): the difference between subsequent values. May be positive (for addition) or negative (for subtraction).

        length (integer): the number of values in the output list."""
        if not isinstance(start, (int, float)):
            raise ValueError("'start' must be an integer or float.")
        if not isinstance(step, (int, float)):
            raise ValueError("'step' must be an integer or float.")
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        output = [start]
        for i in range(1, length):
            denominator = output[i - 1] + step
            if denominator == 0:
                raise ZeroDivisionError("Encountered division by zero in harmonic sequence.")
            output.append(1 / denominator)
        return output
    
    def geometric_seq(self,start,step,length):
        """Return a list of values of the specified length ('length') in which the specified difference between subsequent values ('step') is a ratio (by multiplication).

        Parameters:

        start (integer, float): the first value of the output data range.

        step (integer, float): the difference between subsequent values. May be positive or negative.

        length (integer): the number of values in the output list."""
        if not isinstance(start, (int, float)):
            raise ValueError("'start' must be an integer or float.")
        if not isinstance(step, (int, float)):
            raise ValueError("'step' must be an integer or float.")
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        output = [start]
        for i in range(1, length):
            output.append(output[i - 1] * step)
        return output

    def prime_seq(self, length):
        """Return a list of the first ('length') prime numbers.
        
        Parameters:
        
        length (integer): the number of prime numbers in the output list. """
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        def is_prime(n):
            if n < 2:
                return False
            if n == 2:
                return True
            if n % 2 == 0:
                return False
            for i in range(3, int(n**0.5) + 1, 2):
                if n % i == 0:
                    return False
            return True
        output = []
        num = 2
        while len(output) < length:
            if is_prime(num):
                output.append(num)
            num += 1
        return output

    def fibonacci_seq(self, length):
        """Return a list of the first ('length') Fibonacci numbers.
        
        Parameters:
        
        length (integer): the number of Fibonacci numbers in the output list."""
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        output = [0, 1]
        while len(output) < length:
            output.append(output[-1] + output[-2])
        return output[:length]

    # ============================================================
    # NUMERIC GENERATORS
    # ============================================================

    def rand_ints(self,min_val,max_val,length):
        """Return a list of random integers of the specified length ('length') in the inclusive range specified by the minimum ('min_val') and maximum ('max_val') values.

        Parameters:

        min_val (integer): the minimum value of the output data range.

        max_val (integer): the maximum value of the output data range.

        length (integer): the number of values in the output list."""
        if not isinstance(min_val, int):
            raise ValueError("'min_val' must be an integer.")
        if not isinstance(max_val, int):
            raise ValueError("'max_val' must be an integer.")
        if min_val > max_val:
            raise ValueError("'min_val' must be less than or equal to 'max_val'.")
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        return [random.randint(min_val, max_val) for _ in range(length)]
    
    def rand_floats(self,min_val,max_val,length,decimals = 2):
        """Return a list of random floating-point values of the specified length ('length') in the inclusive range specified by the minimum ('min_val') and maximum ('max_val') values.

        Parameters:

        min_val (integer): the minimum value of the output data range.

        max_val (integer): the maximum value of the output data range.

        length (integer): the number of values in the output list.
        
        decimals (integer): the maximum number of decimal points each value can have."""
        if not isinstance(min_val, (int, float)):
            raise ValueError("'min_val' must be an integer or float.")
        if not isinstance(max_val, (int, float)):
            raise ValueError("'max_val' must be an integer or float.")
        if min_val > max_val:
            raise ValueError("'min_val' must be less than or equal to 'max_val'.")
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")
        return [round(random.uniform(min_val, max_val), decimals) for _ in range(length)]
    
    def rand_points(self, extent, num_points, decimals=2):
        """Return a PointDataset of the specified number of points ('num_points') within the specified extent ('extent').

        Parameters:

        extent (list, tuple): a list or tuple representing the minimum x, minimum y, maximum x, and maximum y, in that order, of the extent to generate points in. 

        num_points (integer): the number of points in the output PointDataset.
        
        decimals (integer): the maximum number of decimal points each coordinate can have.
        """
        if not isinstance(extent, (list, tuple)) or len(extent) != 4:
            raise ValueError("'extent' must be a list or tuple of four values: [x_min, x_max, y_min, y_max].")
        if not all(isinstance(v, (int, float)) for v in extent):
            raise ValueError("All values in 'extent' must be integers or floats.")
        x_min, y_min, x_max, y_max = extent
        if not isinstance(num_points, int) or num_points < 1:
            raise ValueError("'num_points' must be an integer greater than or equal to 1.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")

        points = []
        for i in range(num_points):
            x = round(random.uniform(x_min, x_max), decimals)
            y = round(random.uniform(y_min, y_max), decimals)
            points.append(_Point(i + 1, x, y, []))  # no extra attrs

        headers = ["id", "x", "y"]
        return _PointDataset(points, headers, id_col="id", x_col="x", y_col="y")
    
    def points_in_grid(self, extent, dist, decimals=2):
        """Return a PointDataset within the specified extent ('extent') in a grid pattern. The specified distance between points ('dist') and the extent determine how many points will be in the output PointDataset.
        
        Parameters:

        extent (list, tuple): a list or tuple representing the minimum x, minimum y, maximum x, and maximum y, in that order, of the extent to generate points in. Points are generated from left to right and top to bottom within the extent.

        dist (integer, float): the x and y spacing between points. 
        
        decimals (integer): the maximum number of decimal points each coordinate can have."""
        if not isinstance(extent, (list,tuple)) or len(extent) != 4:
            raise ValueError("'extent' must be a list or tuple of four values: [x_min, x_max, y_min, y_max].")
        if not all(isinstance(v, (int, float)) for v in extent):
            raise ValueError("All values in 'extent' must be integers or floats.")
        x_min, y_min, x_max, y_max = extent
        if not isinstance(dist, (int, float)) or dist <= 0:
            raise ValueError("'dist' must be a positive integer or float.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")

        rows = int((y_max - y_min) / dist) + 1
        cols = int((x_max - x_min) / dist) + 1
        points = []
        count = 1

        for r in range(rows):
            for c in range(cols):
                x = (c * dist) + x_min
                y = y_max - (r * dist)
                if x <= x_max and y >= y_min:
                    points.append(_Point(count, round(x, decimals), round(y, decimals), []))
                    count += 1

        headers = ["id", "x", "y"]
        return _PointDataset(points, headers, id_col="id", x_col="x", y_col="y")

    # ============================================================
    # BOOLEAN, STRINGS, AND COLOURS
    # ============================================================

    def rand_boolean(self,length,numeric = True):
        """Return a list of random Boolean values of the specified length ('length'), which may be either numeric (1 and 0) or not (True and False).

        Parameters:

        length (integer): the number of values in the output data range.

        numeric (Boolean): flag to indicate if the output Boolean values will be represented numerically (1 and 0) or not (True and False)."""
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        if not isinstance(numeric, bool):
            raise ValueError("'numeric' must be a Boolean value.")
        if numeric:
            return [random.randint(0, 1) for _ in range(length)]
        else:
            return [bool(random.randint(0, 1)) for _ in range(length)]

    def rand_strings(self,min_len,max_len,length,chars = "all"):
        """Return a list of the specified length ('length') of random strings between the specified minimum and maximum string length range ('min_len','max_len') using a set of random characters.

        Parameters:

        min_len (integer): the minimum length of each string in the output list.

        max_len (integer): the maximum length of each string in the output list.

        length (integer): the number of strings in the output list.
        
        chars (string, list): a set of characters to build random strings from. The string "all" may be specified to use all possible characters. Alternatively, a list may be specified to use only select characters, including "Upper" for upper-case alphabetic characters, "Lower" for lower-case alphabetic characters, "Numeric" for numeric characters, and "Special" for special characters such as punctuation."""
        import string

        if not isinstance(min_len, int) or min_len < 1:
            raise ValueError("'min_len' must be an integer greater than or equal to 1.")
        if not isinstance(max_len, int) or max_len < min_len:
            raise ValueError("'max_len' must be an integer greater than or equal to 'min_len'.")
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        if not (isinstance(chars, str) or isinstance(chars, list)):
            raise ValueError("'chars' must be either 'all' or a list of character categories.")
        upper = list(string.ascii_uppercase)
        lower = list(string.ascii_lowercase)
        numeric = list(string.digits)
        special = list(string.punctuation)
        use_chars = []
        if chars == "all":
            use_chars = upper + lower + numeric + special
        else:
            if "Upper" in chars:
                use_chars += upper
            if "Lower" in chars:
                use_chars += lower
            if "Numeric" in chars:
                use_chars += numeric
            if "Special" in chars:
                use_chars += special
            if not use_chars:
                raise ValueError("No valid character categories specified in 'chars'.")
        output = []
        for _ in range(length):
            str_len = random.randint(min_len, max_len)
            string_to_add = "".join(random.choice(use_chars) for _ in range(str_len))
            output.append(string_to_add)
        return output

    def rand_strings_user_defined(self,min_len,max_len,length,chars):
        """Return a list of the specified length ('length') of random strings between the specified minimum and maximum string length range ('min_len','max_len') using a set of user-defined characters ('chars').

        Parameters:

        min_len (integer): the minimum length of each string in the output list.

        max_len (integer): the maximum length of each string in the output list.

        length (integer): the number of strings in the output list.
        
        chars (list): a list of characters to build random strings from."""
        if not isinstance(min_len, int) or min_len < 1:
            raise ValueError("'min_len' must be an integer greater than or equal to 1.")
        if not isinstance(max_len, int) or max_len < min_len:
            raise ValueError("'max_len' must be an integer greater than or equal to 'min_len'.")
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        if not isinstance(chars, list) or not chars:
            raise ValueError("'chars' must be a non-empty list of characters.")
        output = []
        for _ in range(length):
            str_len = random.randint(min_len, max_len)
            string_to_add = "".join(random.choice(chars) for _ in range(str_len))
            output.append(string_to_add)
        return output

    def rand_colours(self,length,fmt="rgb"):
        """Return a list of the specified length ('length') of random colours of the specified colour space format ('fmt').
        
        Parameters:
        
        length (integer): the number of values in the output list.
        
        fmt (string): the colour space format of the output colours. May be specified as "rgb" for RGB (Red Green Blue), "cmy" for CMY (Cyan Magenta Yellow), "hex" for hex code, "hsv" for HSV (Hue Saturation Value), and "hsl" for HSL (Hue Saturation Lightness)."""
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        if fmt not in {"rgb", "cmy", "hex", "hsv", "hsl"}:
            raise ValueError("'fmt' must be one of: 'rgb', 'cmy', 'hex', 'hsv', 'hsl'.")
        if fmt == "rgb":
            return [(random.randint(0, 255), random.randint(0, 255), random.randint(0, 255)) for _ in range(length)]
        elif fmt == "cmy":
            return [
                (
                    round(1 - (random.randint(0, 255) / 255), 2),
                    round(1 - (random.randint(0, 255) / 255), 2),
                    round(1 - (random.randint(0, 255) / 255), 2),
                )
                for _ in range(length)
            ]
        elif fmt == "hex":
            return [
                "#%02x%02x%02x" % (random.randint(0, 255), random.randint(0, 255), random.randint(0, 255))
                for _ in range(length)
            ]
        elif fmt == "hsv":
            return [(random.randint(0, 360), round(random.random(), 2), round(random.random(), 2)) for _ in range(length)]
        elif fmt == "hsl":
            return [(random.randint(0, 360), round(random.random(), 2), round(random.random(), 2)) for _ in range(length)]

    # ============================================================
    # CONTINUOUS PROBABILITY DISTRIBUTIONS
    # ============================================================    
    
    def gaussian(self,mean,sd,length,decimals = 2):
        """Return a list of values of the specified length ('length') to fit a Gaussian distribution.

        Parameters:

        mean (integer, float): the mean of the distribution.

        sd (integer, float): the standard deviation of the distribution.

        length (integer): the number of values in the output list.

        decimals (integer): the number of decimal places values will be rounded to."""
        if not isinstance(mean, (int, float)):
            raise ValueError("'mean' must be an integer or float.")
        if not isinstance(sd, (int, float)) or sd <= 0:
            raise ValueError("'sd' must be a positive integer or float.")
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")
        return [round(random.gauss(mean, sd), decimals) for _ in range(length)]

    def triangular(self,min_val,max_val,mode,length,decimals = 2):
        """Return a list of values of the specified length ('length') to fit a triangular distribution (approximately symmetric about the mode) within the specified range.
        
        Parameters:

        min_val (integer): the minimum value of the output data range.

        max_val (integer): the maximum value of the output data range.

        mode (integer, float): the mode of the distribution.
        
        length (integer): the number of values in the output list.

        decimals (integer): the number of decimal places values will be rounded to."""
        if not isinstance(min_val, (int, float)):
            raise ValueError("'min_val' must be an integer or float.")
        if not isinstance(max_val, (int, float)):
            raise ValueError("'max_val' must be an integer or float.")
        if not isinstance(mode, (int, float)):
            raise ValueError("'mode' must be an integer or float.")
        if min_val > max_val:
            raise ValueError("'min_val' must be less than or equal to 'max_val'.")
        if not (min_val <= mode <= max_val):
            raise ValueError("'mode' must be between 'min_val' and 'max_val'.")
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")
        return [round(random.triangular(min_val, max_val, mode), decimals) for _ in range(length)]

    def beta(self,alpha,beta,length,decimals =  2):
        """Return a list of values of the specified length ('length') to fit a beta distribution.
        
        Parameters:

        alpha (integer, float): the alpha value of the output distribution. Must be greater than 0.

        beta (integer, float): the beta value of the output distribution. Must be greater than 0.
        
        length (integer): the number of values in the output list.

        decimals (integer): the number of decimal places values will be rounded to."""
        if not isinstance(alpha, (int, float)) or alpha <= 0:
            raise ValueError("'alpha' must be a positive integer or float.")
        if not isinstance(beta, (int, float)) or beta <= 0:
            raise ValueError("'beta' must be a positive integer or float.")
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")
        return [round(random.betavariate(alpha, beta), decimals) for _ in range(length)]

    def exponential(self,lmbda,length,decimals = 2):
        """Return a list of values of the specified length ('length') to fit an exponential distribution.
        
        Parameters:

        lmbda (integer, float): the lambda value of the distribution. Must not be 0. The lambda value should be the reciprocal of the desired mean.
        
        length (integer): the number of values in the output list.

        decimals (integer): the number of decimal places values will be rounded to."""
        if not isinstance(lmbda, (int, float)) or lmbda == 0:
            raise ValueError("'lmbda' must be a non-zero integer or float.")
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")
        return [round(random.expovariate(lmbda), decimals) for _ in range(length)]
    
    def gamma(self,alpha,beta,length,decimals = 2):
        """Return a list of values of the specified length ('length') to fit a gamma distribution.
        
        Parameters:

        alpha (integer, float): the alpha value of the output distribution. Must be greater than 0.

        beta (integer, float): the beta value of the output distribution. Must be greater than 0.
        
        length (integer): the number of values in the output list.

        decimals (integer): the number of decimal places values will be rounded to."""
        if not isinstance(alpha, (int, float)) or alpha <= 0:
            raise ValueError("'alpha' must be a positive integer or float.")
        if not isinstance(beta, (int, float)) or beta <= 0:
            raise ValueError("'beta' must be a positive integer or float.")
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")
        return [round(random.gammavariate(alpha, beta), decimals) for _ in range(length)]
    
    def lognormal(self,mean,sd,length,decimals = 2):
        """Return a list of values of the specified length ('length') to fit a log normal distribution.

        Parameters:

        mean (integer, float): the mean of the distribution.

        sd (integer, float): the standard deviation of the distribution. Must be greater than 0.

        length (integer): the number of values in the output list.
        
        decimals (integer): the number of decimal places values will be rounded to."""
        if not isinstance(mean, (int, float)):
            raise ValueError("'mean' must be an integer or float.")
        if not isinstance(sd, (int, float)) or sd <= 0:
            raise ValueError("'sd' must be a positive integer or float.")
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")
        return [round(random.lognormvariate(mean, sd), decimals) for _ in range(length)]

    def pareto(self,alpha,length,decimals = 2):
        """Return a list of values of the specified length ('length') to fit a Pareto distribution.
        
        Parameters:

        alpha (integer, float): the alpha value (shape) of the output distribution.
        
        length (integer): the number of values in the output list.

        decimals (integer): the number of decimal places values will be rounded to."""
        if not isinstance(alpha, (int, float)) or alpha <= 0:
            raise ValueError("'alpha' must be a positive integer or float.")
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")
        return [round(random.paretovariate(alpha), decimals) for _ in range(length)]
    
    def weibull(self,alpha,beta,length,decimals = 2):
        """Return a list of values of the specified length ('length') to fit a Weibull distribution.
        
        Parameters:

        alpha (integer, float): the alpha value (scale) of the output distribution.

        beta (integer, float): the beta value (shape) of the output distribution.
        
        length (integer): the number of values in the output list.

        decimals (integer): the number of decimal places values will be rounded to."""
        if not isinstance(alpha, (int, float)) or alpha <= 0:
            raise ValueError("'alpha' must be a positive integer or float.")
        if not isinstance(beta, (int, float)) or beta <= 0:
            raise ValueError("'beta' must be a positive integer or float.")
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")
        return [round(random.weibullvariate(alpha, beta), decimals) for _ in range(length)]

    def uniform(self, min_val, max_val, length, decimals=2):
        """Return a list of values of the specified length ('length') to fit a uniform distribution.
        
        Parameters:
        
        min_val (integer or float): the minimum value of the output range.
        
        max_val (integer or float): the maximum value of the output range.
        
        length (integer): the number of values in the output list.
        
        decimals (integer): the number of decimal places values will be rounded to."""
        if not isinstance(min_val, (int, float)):
            raise ValueError("'min_val' must be an integer or float.")
        if not isinstance(max_val, (int, float)):
            raise ValueError("'max_val' must be an integer or float.")
        if min_val > max_val:
            raise ValueError("'min_val' must be less than or equal to 'max_val'.")
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")
        return [round(random.uniform(min_val, max_val), decimals) for _ in range(length)]

    # ============================================================
    # DISCRETE PROBABILITY DISTRIBUTIONS
    # ============================================================

    def binomial(self, n, p, length):
        """Return a list of values of the specified length ('length') to fit a binomial distribution.
        
        Parameters:
        
        n (integer): the number of trials. Must be greater than 0.
        
        p (float): the probability of success in each trial. Must be between 0 and 1.
        
        length (integer): the number of values in the output list."""
        if not isinstance(n, int) or n <= 0:
            raise ValueError("'n' must be a positive integer.")
        if not isinstance(p, (int, float)) or not (0 <= p <= 1):
            raise ValueError("'p' must be a float between 0 and 1.")
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        output = []
        for _ in range(length):
            successes = sum(1 for _ in range(n) if random.random() < p)
            output.append(successes)
        return output

    def poisson(self, lam, length):
        """Return a list of values of the specified length ('length') to fit a Poisson distribution.
        
        Parameters:
        
        lam (float): the expected number of events (lambda). Must be greater than 0.
        
        length (integer): the number of values in the output list. """
        if not isinstance(lam, (int, float)) or lam <= 0:
            raise ValueError("'lam' must be a positive integer or float.")
        if not isinstance(length, int) or length < 1:
            raise ValueError("'length' must be an integer greater than or equal to 1.")
        output = []
        for _ in range(length):
            L = math.exp(-lam)
            k = 0
            p = 1.0
            while p > L:
                k += 1
                p *= random.random()
            output.append(k - 1)
        return output

class _Math():
    """A set of functions for performing mathematical operations on lists (or nested lists) of numerical values."""
    # ============================================================
    # PRIVATE METHODS AND CLASS ATTRIBUTES AND FUNCTION DISPLAY
    # ============================================================
    def __init__(self):
        self._total_funcs = len([f for f in dir(__class__) if not f.startswith('__')])
    
    def __repr__(self):
        """If the toolbox is printed, display a message."""
        return ".math: a set of functions for performing mathematical operations on lists (or nested lists) of numerical values."

    def display_toolbox_functions(self):
        """Display a list of all available functions within this toolbox."""
        print(f"Number of {__class__.__name__[1:]} functions: {len([f for f in dir(__class__) if not f.startswith('__')])}")
        for f in [f for f in dir(__class__) if not f.startswith("__")]:
            print(f)

    # ============================================================
    # BASIC ARITHMETIC OPERATIONS
    # ============================================================

    def basic_math(self, vals, operand, mode="add", decimals=2):
        """Perform a basic arithmetic operation on a list of values ('vals'), either with a constant number or pairwise element-wise with another list of equal length and return a new list.

        Parameters:

        vals (list): the list of numeric values the operation will be performed on. The input list or lists will not be altered by the operation.

        operand (integer, float, list): either a constant number to apply to each element in 'vals', or a list of equal length for pairwise element-wise operations.

        mode (string): the arithmetic operation to perform. Must be one of 'add', 'subtract', 'multiply', 'divide'.

        decimals (integer): the number of decimal places float values will be rounded to. """
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        if mode not in ("add", "subtract", "multiply", "divide"):
            raise ValueError("'mode' must be one of: 'add', 'subtract', 'multiply', 'divide'.")

        # Constant mode
        if isinstance(operand, (int, float)):
            b_list = [operand] * len(vals)
        # Pairwise mode
        elif isinstance(operand, list):
            if len(vals) != len(operand):
                raise ValueError("Lists must be of equal length for pairwise operations.")
            if not all(isinstance(v, (int, float)) for v in operand):
                raise ValueError("All elements in 'operand' must be numeric.")
            b_list = operand
        else:
            raise TypeError("'operand' must be either a number or a list of numbers.")

        # Perform operation
        output = []
        for a, b in zip(vals, b_list):
            if mode == "add":
                result = a + b
            elif mode == "subtract":
                result = a - b
            elif mode == "multiply":
                result = a * b
            elif mode == "divide":
                if b == 0:
                    raise ZeroDivisionError("Division by zero encountered.")
                result = a / b
            output.append(round(result, decimals))

        return output

    def average(self,vals,decimals = 2):
        """Calculate the mean of the input list ('vals') and return a single value.
        
        Parameters:
        
        vals (list): the list of values (or nested list) the operation will be performed on. The input list or lists will not be altered by the operation.
        
        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        return round(sum(vals)/len(vals),decimals)

    def sum_squares(self,vals,decimals = 2):
        """Calculate and return the sum of squared deviations from the mean for an input variable.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        mean = sum(vals) / len(vals)
        sum_sqr = sum((val - mean) ** 2 for val in vals)

        return round(sum_sqr, decimals)
    
    def diff_squares(self,x_var,y_var,decimals = 2):
        """Calculate and return the difference of squared deviations between two input variables.

        Parameters:

        x_var (list): the list of values representing the x variable the operation will be performed on. The input list will not be altered by the operation.

        y_var (list): the list of values representing the y variable the operation will be performed on. The input list will not be altered by the operation.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(x_var, list) or not x_var:
            raise TypeError("'x_var' must be a non-empty list.")
        if not isinstance(y_var, list) or not y_var:
            raise TypeError("'y_var' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in x_var):
            raise ValueError("All elements in 'x_var' must be numeric.")
        if not all(isinstance(v, (int, float)) for v in y_var):
            raise ValueError("All elements in 'y_var' must be numeric.")
        if len(x_var) != len(y_var):
            raise ValueError("'x_var' and 'y_var' must be equal length.")
        xmean = sum(x_var) / len(x_var)
        ymean = sum(y_var) / len(y_var)

        sum_sqr_x = sum((val - xmean) ** 2 for val in x_var)
        sum_sqr_y = sum((val - ymean) ** 2 for val in y_var)

        return round(sum_sqr_x - sum_sqr_y, decimals)
    
    def accumulate(self,vals,decimals = 2):
        """For every value in an input list ('vals'), calculate the cumulative sum and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        v_sum = 0
        output = []
        append = output.append
        for val in vals:
            v_sum += val
            append(round(v_sum,decimals))
        return output
    
    def power(self,vals,power,decimals = 2):
        """Raise each value in an input list ('vals') to a specified power ('power') and return a new list.
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        power (integer, float): the exponent to raise each element in the list by.
        
        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        if not isinstance(power, (int,float)):
            raise ValueError("'power' must be numeric.")
        return [round(val**power,decimals) for val in vals]

    def root(self,vals,root,decimals = 2):
        """Calculate the nth root ('root') of each value in an input list ('vals') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        root (integer): the nth root to take of each value.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        if not isinstance(root, (int, float)):
            raise TypeError("'root' must be numeric.")
        output = []
        for val in vals:
            if val < 0 and root % 2 == 0:
                raise ValueError(f"Cannot take an even root of a negative value: {val}")
            output.append(round(val ** (1 / root), decimals))

        return output
 
    def ln(self,vals,decimals = 2):
        """Calculate the natural logarithm of each value in an input list ('vals') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        if any(v <= 0 for v in vals):
            raise ValueError("All elements in 'vals' must be positive to compute natural logarithm.")
        return [round(math.log(val),decimals) for val in vals]
    
    def log2(self,vals,decimals = 2):
        """Calculate the base-2 logarithm of each value in an input list ('vals') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        if any(v <= 0 for v in vals):
            raise ValueError("All elements in 'vals' must be positive to compute base-2 logarithm.")
        return [round(math.log2(val),decimals) for val in vals]
    
    def log10(self,vals,decimals = 2):
        """Calculate the base-10 logarithm of each value in an input list ('vals') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        if any(v <= 0 for v in vals):
            raise ValueError("All elements in 'vals' must be positive to compute base-10 logarithm.")
        return [round(math.log10(val),decimals) for val in vals]

    def reciprocal(self,vals,decimals = 2):
        """Calculate the reciprocal (1/value) of each value in an input list ('vals') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        if any(v == 0 for v in vals):
            raise ValueError("All elements in 'vals' must be non-zero to compute reciprocals.")
        return [round(1/val,decimals) for val in vals]

    def absolute(self,vals,decimals = 2):
        """Calculate the absolute value of each value in an input list ('vals') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        return [round(abs(val),decimals) for val in vals]
    
    def negate(self,vals,decimals = 2):
        """Negate each value in an input list ('vals') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        return [round(-val,decimals) for val in vals]

    # ============================================================
    # ROUNDING AND TRUNCATION
    # ============================================================

    def rounding(self,vals,decimals = 2):
        """Round each value in an input list ('vals') to a specified number of decimal places ('decimals') and return a new list.
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        return [round(val,decimals) for val in vals]

    def ceiling(self,vals):
        """Calculate the ceiling (the smallest integer greater than or equal to) of each value in an input list ('vals') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        return [math.ceil(val) for val in vals]
        
    def floor(self,vals):
        """Calculate the floor (the largest integer less than or equal to) of each value in an input list ('vals') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        return [math.floor(val) for val in vals]

    def trunc(self, vals):
        """Truncate each value in an input list ('vals') toward zero and return a new list.

        Parameters:

        vals (list): the list of numeric values the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        return [math.trunc(val) for val in vals]

    # ============================================================
    # TRIGONOMETRY AND ANGLE CONVERSIONS
    # ============================================================

    def to_radians(self,vals,decimals = 2):
        """Convert each value in an input list ('vals') from degrees to radians and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        return [round(math.radians(val), decimals) for val in vals]

    def to_degrees(self,vals,decimals = 2):
        """Convert each value in an input list ('vals') from radians to degrees and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        return [round(math.degrees(val), decimals) for val in vals]
    
    def trig_ops(self, vals, operation, decimals=2):
        """Apply a trigonometric operation to each value in the input list ('vals') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        operation (string): the trigonometric operation to apply. Possible values include "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh", "tanh", "asinh", "acosh", "atanh".

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")

        valid_ops = {
            "sin": math.sin,
            "cos": math.cos,
            "tan": math.tan,
            "asin": math.asin,
            "acos": math.acos,
            "atan": math.atan,
            "sinh": math.sinh,
            "cosh": math.cosh,
            "tanh": math.tanh,
            "asinh": math.asinh,
            "acosh": math.acosh,
            "atanh": math.atanh
        }

        if operation not in valid_ops:
            raise ValueError("Unsupported trigonometric operation.")
        
        if operation in {"asin", "acos"} and not all(-1 <= v <= 1 for v in vals):
            raise ValueError(f"All values must be between -1 and 1 for {operation}.")
        if operation == "acosh" and not all(v >= 1 for v in vals):
            raise ValueError("All values must be greater than or equal to 1 for acosh.")
        if operation == "atanh" and not all(-1 < v < 1 for v in vals):
            raise ValueError("All values must be strictly between -1 and 1 for atanh.")


        func = valid_ops[operation]
        return [round(func(v), decimals) for v in vals]

    # ============================================================
    # RESCALING AND NORMALIZATION
    # ============================================================

    def linear_rescale(self,vals,new_min,new_max, inverted = False, decimals = 2):
        """Rescale each value in an input list ('vals') to fit within a specified range ('new_min','new_max') and return a new list. The rescale function performs a linear stretch.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        new_min (integer, float): the new minimum value of the data after rescaling.

        new_max (integer, float): the new maximum value of the data after rescaling.

        inverted (Boolean): flag to indicate if the rescaled values will be inverted so that lower values in the input data are rescaled to higher values in the output data and vice versa.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        if not isinstance(new_min, (int, float)) or not isinstance(new_max, (int, float)):
            raise TypeError("'new_min' and 'new_max' must be numeric.")
        if new_min == new_max:
            raise ValueError("'new_min' and 'new_max' must be distinct values.")
        output = []
        max_v = max(vals)
        min_v = min(vals)
        for val in vals:
            if not inverted:
                op = (val - min_v)/(max_v - min_v)
                val = ((new_max - new_min) * op) + new_min
            else:
                op = 1 - ((val - min_v)/(max_v - min_v)) 
                val = ((new_max - new_min) * op) + new_min
            output.append(round(val,decimals))
        return output

    def linear_rescale_categorical(self,vals,variables,new_min,new_max,decimals = 2,inverted = False):
        """Rescale each categorical value (integer or string) in an input list ('vals') to fit within a specified numerical range ('new_min','new_max') and return a new list. The rescale function performs a linear stretch.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        variables (list): the categorical variables to rescale, in the order they will be rescaled (i.e. the first variable will be assigned the new minimum, and the last variable will be assigned the new maximum). Must contain any/all possible values in the input list.

        new_min (integer, float): the new minimum value of the data after rescaling.

        new_max (integer, float): the new maximum value of the data after rescaling.

        inverted (Boolean): flag to indicate if the rescaled values will be inverted so that lower values in the input data are rescaled to higher values in the output data and vice versa.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not isinstance(variables, list) or not variables:
            raise TypeError("'variables' must be a non-empty list.")
        if len(set(variables)) != len(variables):
            raise ValueError("'variables' must contain unique values.")
        if not all(v in variables for v in vals):
            raise ValueError("All elements in 'vals' must be present in 'variables'.")
        if not isinstance(new_min, (int, float)) or not isinstance(new_max, (int, float)):
            raise TypeError("'new_min' and 'new_max' must be numeric.")
        if new_min == new_max:
            raise ValueError("'new_min' and 'new_max' must be distinct values.")
        if len(variables) == 1:
            raise ValueError("'variables' must contain at least two distinct values for rescaling.")

        min_t, max_t = 0, len(variables) - 1
        output = []
        for val in vals:
            idx = variables.index(val)
            if not inverted:
                op = (idx - min_t) / (max_t - min_t)
            else:
                op = 1 - ((idx - min_t) / (max_t - min_t))
            rescaled = ((new_max - new_min) * op) + new_min
            output.append(round(rescaled, decimals))

        return output

    def normalize(self, vals, decimals=2):
        """Normalize the values in the specified list ('vals') to a range between 0 and 1 using min-max scaling and return a new list. Missing or invalid values will be ignored. Normalized values will be returned as floats.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")

        # Filter numeric values
        numeric_vals = [v for v in vals if isinstance(v, (int, float))]
        if not numeric_vals:
            raise ValueError("List has no valid numeric values to normalize.")

        min_val, max_val = min(numeric_vals), max(numeric_vals)
        if min_val == max_val:
            raise ValueError("Cannot normalize when all values are identical.")

        output = []
        for v in vals:
            if isinstance(v, (int, float)):
                norm = (v - min_val) / (max_val - min_val)
                output.append(round(norm, decimals))
            else:
                output.append(v)

        return output

    def min_max_clip(self,vals,min_val = None,max_val = None,decimals = 2):
        """Clip each value in an input list ('vals') to fit within the specified minimum ('min_val') and/or maximum ('max_val') value range and return a new list. Values greater than the maximum or less than the minimum will be set to the specified maximum or specified minimum, respectively. 
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        min_val (integer, float): the minimum value from the input list that will be included in the output list. Input values smaller than the minimum will be set to the minimum. May be left unspecified for a one-tail clip.

        max_val (integer, float): the maximum value from the input list that will be included in the output list. Input values larger than the maximum will be set to the maximum. May be left unspecified for a one-tail clip.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        if min_val is not None and not isinstance(min_val, (int, float)):
            raise TypeError("'min_val' must be numeric if specified.")
        if max_val is not None and not isinstance(max_val, (int, float)):
            raise TypeError("'max_val' must be numeric if specified.")
        if min_val is not None and max_val is not None and min_val > max_val:
            raise ValueError("'min_val' must be less than or equal to 'max_val'.")

        output = []
        for val in vals:
            clipped = val
            if min_val is not None and clipped < min_val:
                clipped = min_val
            if max_val is not None and clipped > max_val:
                clipped = max_val
            output.append(round(clipped, decimals))

        return output

    def percent_clip(self,vals,percent,tail="Both",decimals = 2):
        """Clip each value in an input list ('vals') to fit within the specified percentage ('percent') of either/both tails of the distribution. Values greater than or less than the percentage clip will be set to the maximum or minimum value, respectively.
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        percent (integer): the percentage of values to clip from either/both tails of the distribution.

        tail (string): which tail of the distribution to apply the clip to. A value of "Upper" or "Lower" may be specified for an upper or lower one-tail clip respectively. A value of "Both" will clip values at both ends of the distribution.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        if not isinstance(percent, (int, float)):
            raise TypeError("'percent' must be numeric.")
        if percent < 0 or percent > 100:
            raise ValueError("'percent' must be between 0 and 100.")
        if tail not in {"Upper", "Lower", "Both"}:
            raise ValueError("'tail' must be one of 'Upper', 'Lower', or 'Both'.")

        v_max = max(vals)
        v_min = min(vals)
        val_range = v_max - v_min
        if val_range == 0:
            raise ValueError("Cannot clip when all input values are identical.")

        to_remove = val_range * (percent / 100)
        output = []

        for val in vals:
            if tail == "Upper":
                clipped = val if val <= v_max - to_remove else v_max - to_remove
            elif tail == "Lower":
                clipped = val if val >= v_min + to_remove else v_min + to_remove
            else:  # Both
                if val < v_min + to_remove:
                    clipped = v_min + to_remove
                elif val > v_max - to_remove:
                    clipped = v_max - to_remove
                else:
                    clipped = val
            output.append(round(clipped, decimals))

        return output

    # ============================================================
    # INTERPOLATION AND FILLING
    # ============================================================

    def interpolate_single_value(self, vals, fill_val, to_fill=""):
        """Replace each missing or invalid value ('to_fill') in a list ('vals') with a single specified value ('fill_val') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        fill_val (integer, float, string): the new value that will replace missing or invalid values.

        to_fill (integer, float, string, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")

        # Normalize to_fill into a list
        if not isinstance(to_fill, list):
            to_fill = [to_fill]

        output = []
        for v in vals:
            if v in to_fill:
                output.append(fill_val)
            else:
                output.append(v)

        return output
    
    def interpolate_from_list(self, vals, fill_list, to_fill=""):
        """Replace each missing or invalid value ('to_fill') in a list ('vals') with the corresponding value from another list ('fill_list') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        fill_list (list): the list of values that will be used to fill missing values. Must be the same length as 'vals'.

        to_fill (integer, float, string, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value. """
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not isinstance(fill_list, list) or not fill_list:
            raise TypeError("'fill_list' must be a non-empty list.")
        if len(fill_list) != len(vals):
            raise ValueError("Length of 'fill_list' must be equal to length of 'vals'.")

        # Normalize to_fill into a list
        if not isinstance(to_fill, list):
            to_fill = [to_fill]

        output = []
        for ix, v in enumerate(vals):
            if v in to_fill and fill_list[ix] not in to_fill:
                output.append(fill_list[ix])
            else:
                output.append(v)

        return output
    
    def interpolate_stat(self, vals, stat, to_fill="", decimals=2):
        """Interpolate a new value for each missing or invalid value ('to_fill') in the specified list ('vals') using the specified statistic ('stat') of all valid values within that list and return a new list. Missing or invalid values will be interpolated as floats.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        stat (string): the specific statistic to return. Possible statistics are "Count", "Unique Values", "First", "Last", "Sum", "Max", "Min", "Range", "Mean", "Median", "Mode", "IQ1", "IQ3", "IQR", "Standard Deviation", "Standard Error", "Variance", "Coefficient of Variation", "Skewness", and "Kurtosis".

        to_fill (integer, float, string, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value.

        decimals (integer): the number of decimal places float values will be rounded to. """
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")

        valid_stats = {
            "Count", "Unique Values", "First", "Last", "Sum", "Max", "Min", "Range", "Mean",
            "Median", "Mode", "IQ1", "IQ3", "IQR", "Standard Deviation",
            "Standard Error", "Variance", "Coefficient of Variation",
            "Skewness", "Kurtosis"
        }

        if stat not in valid_stats:
            raise ValueError(f"Unknown statistic '{stat}'.")

        if not isinstance(to_fill, list):
            to_fill = [to_fill]

        # Filter numeric values
        numeric_vals = [v for v in vals if isinstance(v, (int, float))]
        if not numeric_vals:
            raise ValueError("List has no valid numeric values to compute statistic from.")

        # Compute statistic using your my helper
        stat_val = _Stats().descriptive_stat(numeric_vals, stat, decimals)

        output = []
        for v in vals:
            if v in to_fill:
                output.append(stat_val)
            else:
                output.append(v)

        return output
    
    def interpolate_forward(self, vals, to_fill=""):
        """Interpolate a new value for each missing or invalid value ('to_fill') in the specified list ('vals') using the valid value immediately preceding each missing or invalid value and return a new list. This method will not be able to fill missing or invalid values in the first position.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        to_fill (int, float, str, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")

        if not isinstance(to_fill, list):
            to_fill = [to_fill]

        output = vals[:]
        for i in range(1, len(output)):
            if output[i] in to_fill:
                prev_val = output[i - 1]
                if prev_val not in to_fill:
                    output[i] = prev_val

        return output

    def interpolate_backward(self, vals, to_fill=""):
        """Interpolate a new value for each missing or invalid value ('to_fill') in the specified list ('vals') using the valid value immediately following each missing or invalid value and return a new list. This method will not be able to fill missing or invalid values in the last position.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        to_fill (int, float, str, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")

        if not isinstance(to_fill, list):
            to_fill = [to_fill]

        output = vals[:]
        for i in range(len(output) - 1):
            if output[i] in to_fill:
                next_val = output[i + 1]
                if next_val not in to_fill:
                    output[i] = next_val

        return output
    
    def interpolate_linear(self, vals, to_fill="", direction="forward"):
        """Interpolate a new value for each missing or invalid value ('to_fill') in the specified list ('vals') using a linear interpolation technique and return a new list. This method will not be able to fill missing or invalid values at either the beginning or end of the list, depending on the direction of interpolation.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        to_fill (int, float, str, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value.

        direction (str): the direction of interpolation applied after the linear interpolation. May be specified as "forward" for a forward interpolation (missing or invalid values at the beginning of the list will not be interpolated) or "backward" for a backward interpolation (missing or invalid values at the end of the list will not be interpolated)."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")

        if not isinstance(to_fill, list):
            to_fill = [to_fill]

        rows = len(vals)
        x_col = list(range(1, rows + 1))
        output = vals[:]

        # Perform interpolation
        for i in range(rows):
            if output[i] in to_fill and 0 < i < rows - 1:
                # Find next valid value
                y2, x2 = None, None
                for j in range(i + 1, rows):
                    try:
                        y2 = float(output[j])
                        x2 = x_col[j]
                        break
                    except (ValueError, TypeError):
                        continue

                if y2 is not None:
                    try:
                        x, x1 = x_col[i], x_col[i - 1]
                        y1 = float(output[i - 1])
                        y = y1 + ((x - x1) / (x2 - x1)) * (y2 - y1)
                        output[i] = y
                    except (ValueError, TypeError, ZeroDivisionError):
                        pass

        # Apply forward/backward fill for edge cases
        if direction == "forward":
            output = self.interpolate_forward(output, to_fill)
        elif direction == "backward":
            output = self.interpolate_backward(output, to_fill)

        return output
    
    def interpolate_mean_distance(self, vals, distance, to_fill="", decimals=2):
        """Interpolate a new value for each missing or invalid value ('to_fill') in the specified list ('vals') using the average of valid values before and after each missing or invalid value and return a new list. Missing or invalid values will be interpolated as floats. This technique is best suited to ordered data where the order is meaningful (e.g. data ordered by date or time).

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        distance (integer): the number of values before and after each missing or invalid value to use for interpolation. Only valid (numeric) values within this range will be used for interpolation.

        to_fill (integer, float, string, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not isinstance(distance, int) or distance <= 0:
            raise ValueError("Distance value must be a positive integer.")

        if not isinstance(to_fill, list):
            to_fill = [to_fill]

        output = vals[:]
        rows = len(output)

        for i in range(rows):
            if output[i] in to_fill:
                vals_for_mean = []
                for d in range(1, distance + 1):
                    # Look forward
                    if i + d < rows:
                        val = output[i + d]
                        if val not in to_fill:
                            try:
                                vals_for_mean.append(float(val))
                            except (ValueError, TypeError):
                                pass
                    # Look backward
                    if i - d >= 0:
                        val = output[i - d]
                        if val not in to_fill:
                            try:
                                vals_for_mean.append(float(val))
                            except (ValueError, TypeError):
                                pass

                if vals_for_mean:
                    output[i] = round(sum(vals_for_mean) / len(vals_for_mean), decimals)

        return output
    
    def interpolate_rolling_average(self, vals, window, to_fill="", decimals=2):
        """Interpolate a new value for each missing or invalid value ('to_fill') in the specified list ('vals') using the average of valid values within a rolling window of previous values and return a new list. Missing or invalid values will be interpolated as floats. This technique is best suited to ordered data where the order is meaningful (e.g. data ordered by date or time).

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        window (integer): the number of previous values to use for interpolation. Only valid (numeric) values within this range will be used for interpolation.

        to_fill (integer, float, string, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not isinstance(window, int) or window <= 0:
            raise ValueError("Window size must be a positive integer.")

        if not isinstance(to_fill, list):
            to_fill = [to_fill]

        output = vals[:]
        rows = len(output)

        for i in range(rows):
            if output[i] in to_fill:
                vals_for_mean = []
                for j in range(max(0, i - window), i):  # look only backward
                    val = output[j]
                    if val not in to_fill:
                        try:
                            vals_for_mean.append(float(val))
                        except (ValueError, TypeError):
                            pass
                if vals_for_mean:
                    output[i] = round(sum(vals_for_mean) / len(vals_for_mean), decimals)

        return output

    def interpolate_inverse_distance_weighted(self, vals, distance, power=2, to_fill="", decimals=2):
        """Interpolate a new value for each missing or invalid value ('to_fill') in the specified list ('vals') using an inverse distance weighted (IDW) technique and return a new list. Missing or invalid values will be interpolated as floats. This technique is best suited to ordered data where the order is meaningful (e.g. data ordered by date or time).

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        distance (integer): the number of values before and after each missing or invalid value to use for interpolation. Only valid (numeric) values within this range will be used for interpolation.

        power (integer): the power of the inverse distance weighting function. Higher values will give greater weight to nearer values i.e. nearer values will contribute more strongly to the interpolated value.

        to_fill (integer, float, string, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not isinstance(distance, int) or distance <= 0:
            raise ValueError("Distance value must be a positive integer.")
        if not isinstance(power, int) or power <= 0:
            raise ValueError("Power value must be a positive integer.")

        if not isinstance(to_fill, list):
            to_fill = [to_fill]

        output = vals[:]
        rows = len(output)

        for i in range(rows):
            if output[i] in to_fill:
                numerators = []
                denominators = []
                for d in range(1, distance + 1):
                    # Look forward
                    if i + d < rows:
                        val = output[i + d]
                        if val not in to_fill:
                            try:
                                numerators.append(float(val) / (d ** power))
                                denominators.append(1 / (d ** power))
                            except (ValueError, TypeError):
                                pass
                    # Look backward
                    if i - d >= 0:
                        val = output[i - d]
                        if val not in to_fill:
                            try:
                                numerators.append(float(val) / (d ** power))
                                denominators.append(1 / (d ** power))
                            except (ValueError, TypeError):
                                pass

                if numerators and denominators:
                    output[i] = round(sum(numerators) / sum(denominators), decimals)

        return output
    
    def interpolate_regression(self, vals, x_vals, to_fill="", decimals=2):
        """Interpolate a new value for each missing or invalid value ('to_fill') in the specified list ('vals')
        by performing a linear regression between the values of the specified list ('vals') and a second list ('x_vals') and return a new list. Missing or invalid values will be interpolated as floats.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        x_vals (list): the list of values to perform a linear regression with for interpolation. Only valid (numeric) values within this list will be used for interpolation. Must be the same length as 'vals'.

        to_fill (integer, float, string, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not isinstance(x_vals, list):
            raise TypeError("'x_vals' must be a list.")
        if len(x_vals) != len(vals):
            raise ValueError("Length of 'x_vals' must be equal to length of 'vals'.")

        if not isinstance(to_fill, list):
            to_fill = [to_fill]

        y_vals, x_clean = [], []
        for i in range(len(vals)):
            val = vals[i]
            if val not in to_fill:
                try:
                    y_vals.append(float(val))
                    x_clean.append(float(x_vals[i]))
                except (ValueError, TypeError):
                    pass

        if not y_vals or not x_clean:
            raise ValueError("No valid numeric values available for regression.")

        # Calculate regression coefficients (intercept, slope)
        model = _Stats().linear_regression(x_vals, y_vals)
        intercept, slope = model["intercept"],model["slope"]

        output = vals[:]
        for i in range(len(output)):
            if output[i] in to_fill and x_vals[i] not in to_fill:
                try:
                    output[i] = round((float(x_vals[i]) * slope) + intercept, decimals)
                except (ValueError, TypeError):
                    pass

        return output

    # ============================================================
    # THRESHOLDING AND CONDITIONAL OPERATIONS
    # ============================================================

    def condition_op(self, vals, list_or_val, op="lt", numeric=True):
        """Evaluate a conditional operation element-wise on a list of numeric values and return a new list of Boolean values.

        Parameters:

        vals (list): the list of numeric  values the operation will be performed on. The input list will not be altered by the operation.

        list_or_val (int, float, list): either a single threshold value or a list of values of the same length as 'vals'.

        op (string): the conditional operation to apply. May be specified as 'lt' (<), 'gt' (>), 'le' (<=), 'ge' (>=), 'eq' (==), or 'ne' (!=).

        numeric (bool): flag to indicate if the output values will be represented numerically (1 and 0) or as Boolean (True and False). Default is True."""

        # Validation
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All values in 'vals' must be numeric.")
        if isinstance(list_or_val, list):
            if len(list_or_val) != len(vals):
                raise ValueError("If 'list_or_val' is a list, it must match the length of 'vals'.")
            targets = list_or_val
        elif isinstance(list_or_val, (int, float)):
            targets = [list_or_val] * len(vals)
        else:
            raise ValueError("'list_or_val' must be a numeric value or a list of numeric values.")

        # Operator map
        ops = {
            "lt": lambda a, b: a < b,
            "gt": lambda a, b: a > b,
            "le": lambda a, b: a <= b,
            "ge": lambda a, b: a >= b,
            "eq": lambda a, b: a == b,
            "ne": lambda a, b: a != b,
        }
        if op not in ops:
            raise ValueError(f"Unknown operation '{op}'. Must be one of {list(ops.keys())}.")

        # Apply operation
        result = [ops[op](a, b) for a, b in zip(vals, targets)]

        # Return numeric or Boolean
        return [int(r) if numeric else r for r in result]

    def remove_threshold_condition(self,vals,thresh,condition):
        """Remove values from an input list ('vals') based on an equality threshold condition and return a new list.
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        thresh (integer, float): the threshold to evaluate each value in the list against.
        
        condition (string): the condition placed upon each value in the list to determine if that value will be removed. May be specified as "<" for less than, "<=" for less than equal to, ">" for greater than, ">=" for greater than equal to, "==" for equal to, or "!=" for not equal to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        if not isinstance(thresh, (int, float)):
            raise TypeError("'thresh' must be numeric.")

        ops = {'<':operator.lt, '<=':operator.le, '>': operator.gt, '>=': operator.ge, '==': operator.eq, '!=': operator.ne}
        if condition not in ops:
            raise ValueError("Unsupported condition. Must be one of '<', '<=', '>', '>=', '==', '!='.")

        return [v for v in vals if ops[condition](v,thresh)]
             
    def min_above_threshold(self,vals,thresh):
        """Calculate and return the minimum value of an input list ('vals') greater than a specified threshold ('thresh').
        
        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        thresh (integer, float): the threshold for determining the minimum value of the input list. Only values greater than the threshold will be used to calculate the minimum."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        if not isinstance(thresh, (int, float)):
            raise TypeError("'thresh' must be numeric.")
        return min(self.remove_threshold_condition(vals,thresh,"<="))

    def max_above_threshold(self,vals,thresh):
        """Calculate and return the maximum value of an input list ('vals') greater than a specified threshold ('thresh').
        
        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        thresh (integer, float): the threshold for determining the maximum value of the input list. Only values greater than the threshold will be used to calculate the maximum."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        if not isinstance(thresh, (int, float)):
            raise TypeError("'thresh' must be numeric.")
        return max(self.remove_threshold_condition(vals,thresh,"<="))

    def min_below_threshold(self,vals,thresh):
        """Calculate and return the minimum value of an input list ('vals') less than a specified threshold ('thresh').
        
        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        thresh (integer, float): the threshold for determining the minimum value of the input list. Only values less than the threshold will be used to calculate the minimum."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        if not isinstance(thresh, (int, float)):
            raise TypeError("'thresh' must be numeric.")
        return min(self.remove_threshold_condition(vals,thresh,">="))

    def max_below_threshold(self,vals,thresh):
        """Calculate and return the maximum value of an input list ('vals') less than a specified threshold ('thresh').
        
        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        thresh (integer, float): the threshold for determining the maximum value of the input list. Only values less than the threshold will be used to calculate the maximum."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        if not isinstance(thresh, (int, float)):
            raise TypeError("'thresh' must be numeric.")
        return max(self.remove_threshold_condition(vals,thresh,">="))

    # ============================================================
    # CLEANING
    # ============================================================

    def modify_consecutive_duplicates(self,vals,mod_val,method="add"):
        """Modify consecutive duplicate values in an input list ('vals') by a specified value ('mod_val') and return a new list.
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        mod_val (integer, float): the value to modify consecutive duplicates by.
        
        method (string): the method of modifying consecutive values. May be specified as 'add' to add the modifying value to consecutive duplicates, or 'subtract' to subtract the modifying value from consecutive duplicates."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        if not isinstance(mod_val, (int, float)):
            raise TypeError("'mod_val' must be numeric.")
        if method not in {"add", "subtract"}:
            raise ValueError("'method' must be either 'add' or 'subtract'.")

        data = vals[:]
        for i in range(1, len(data)):
            if data[i] == data[i - 1]:
                if method == "add":
                    data[i] = data[i] + mod_val
                else:  # subtract
                    data[i] = data[i] - mod_val

        return data

    def remove_non_numeric(self,vals):
        """Remove non-numeric values from an input list ('vals') and return a new list.
        
        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        return [v for v in vals if isinstance(v, (int, float))]

    # ============================================================
    # SEQUENCE AND PARITY OPERATIONS
    # ============================================================

    def next_odd(self,vals,direction="Up"):
        """Convert each even value in an input list ('vals') to the next odd number and return a new list. Odd numbers in the input list will not be modified in the output list.

        Parameters:

        vals (list): the list of values the operation will be performed on. Input values should only be integer values. The input list will not be altered by the operation.

        direction (string): the direction of conversion. A value of "Up" indicates that even numbers will be rounded up to the nearest odd number, while a value of "Down" indicates that even numbers will be rounded down to the nearest odd number."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, int) for v in vals):
            raise ValueError("All elements in 'vals' must be integers.")
        if direction not in {"Up", "Down"}:
            raise ValueError("'direction' must be either 'Up' or 'Down'.")

        output = []
        for val in vals:
            if val % 2 == 0:  # even
                if direction == "Up":
                    output.append(val + 1)
                else:  # Down
                    output.append(val - 1)
            else:  # odd
                output.append(val)

        return output

    def next_even(self,vals,direction="Up"):
        """Convert each odd value in an input list ('vals') to the next even number and return a new list. Even numbers in the given list will not be modified in the output list.

        Parameters:

        vals (list): the list of values the operation will be performed on. Input values should only be integer values. The input list will not be altered by the operation.

        direction (string): the direction of conversion. A value of "Up" indicates that odd numbers will be rounded up to the nearest even number, while a value of "Down" indicates that odd numbers will be rounded down to the nearest even number."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, int) for v in vals):
            raise ValueError("All elements in 'vals' must be integers.")
        if direction not in {"Up", "Down"}:
            raise ValueError("'direction' must be either 'Up' or 'Down'.")

        output = []
        for val in vals:
            if val % 2 == 1:  # odd
                if direction == "Up":
                    output.append(val + 1)
                else:  # Down
                    output.append(val - 1)
            else:  # even
                output.append(val)

        return output

    def is_even_odd(self, vals, mode="even", numeric=True):
        """Evaluate if each numeric value in the input list ('vals') is even or odd and return a list of Boolean values.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        mode (string): specify 'even' or 'odd' to determine which check to perform.

        numeric (bool): flag to indicate if the output values will be represented numerically (1 and 0) or as Boolean (True and False). Default is True."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if mode not in {"even", "odd"}:
            raise ValueError("Mode must be 'even' or 'odd'.")
        if not all(isinstance(v, int) for v in vals):
            raise TypeError("All elements in 'vals' must be integers.")

        if mode == "even":
            result = [val % 2 == 0 for val in vals]
        else:  # odd
            result = [val % 2 == 1 for val in vals]
        return [int(r) if numeric else r for r in result]
 
class _Stats():
    """A set of functions for performing statistical operations on lists of numerical values."""
    # ============================================================
    # PRIVATE METHODS AND CLASS ATTRIBUTES AND FUNCTION DISPLAY
    # ============================================================
    def __init__(self):
        self._total_funcs = len([f for f in dir(__class__) if not f.startswith('_')])
    
    _T_CRITICAL_VALUES = {
        # Degrees of Freedom (df):
        1: { 
            0.30: {'two_tail': 1.895, 'one_tail': 1.376}, 0.25: {'two_tail': 2.414, 'one_tail': 1.558}, 
            0.20: {'two_tail': 3.078, 'one_tail': 1.963}, 0.15: {'two_tail': 4.160, 'one_tail': 2.571},
            0.10: {'two_tail': 6.314, 'one_tail': 3.143}, 0.05: {'two_tail': 12.706, 'one_tail': 6.314}, 
            0.01: {'two_tail': 63.657, 'one_tail': 31.821}
        },
        2: { 
            0.30: {'two_tail': 1.386, 'one_tail': 1.061}, 0.25: {'two_tail': 1.604, 'one_tail': 1.190}, 
            0.20: {'two_tail': 1.886, 'one_tail': 1.386}, 0.15: {'two_tail': 2.353, 'one_tail': 1.706},
            0.10: {'two_tail': 2.920, 'one_tail': 1.886}, 0.05: {'two_tail': 4.303, 'one_tail': 2.920}, 
            0.01: {'two_tail': 9.925, 'one_tail': 6.965}
        },
        3: {
            0.30: {'two_tail': 1.250, 'one_tail': 0.978}, 0.25: {'two_tail': 1.423, 'one_tail': 1.108}, 
            0.20: {'two_tail': 1.638, 'one_tail': 1.250}, 0.15: {'two_tail': 1.996, 'one_tail': 1.533},
            0.10: {'two_tail': 2.353, 'one_tail': 1.943}, 0.05: {'two_tail': 3.182, 'one_tail': 2.353}, 
            0.01: {'two_tail': 5.841, 'one_tail': 4.541}
        },
        4: {
            0.30: {'two_tail': 1.190, 'one_tail': 0.941}, 0.25: {'two_tail': 1.344, 'one_tail': 1.054}, 
            0.20: {'two_tail': 1.533, 'one_tail': 1.190}, 0.15: {'two_tail': 1.826, 'one_tail': 1.439},
            0.10: {'two_tail': 2.132, 'one_tail': 1.638}, 0.05: {'two_tail': 2.776, 'one_tail': 2.132}, 
            0.01: {'two_tail': 4.604, 'one_tail': 3.747}
        },
        5: {
            0.30: {'two_tail': 1.156, 'one_tail': 0.920}, 0.25: {'two_tail': 1.301, 'one_tail': 1.025}, 
            0.20: {'two_tail': 1.476, 'one_tail': 1.156}, 0.15: {'two_tail': 1.746, 'one_tail': 1.376},
            0.10: {'two_tail': 2.015, 'one_tail': 1.476}, 0.05: {'two_tail': 2.571, 'one_tail': 2.015}, 
            0.01: {'two_tail': 4.032, 'one_tail': 3.365}
        },
        6: {
            0.30: {'two_tail': 1.134, 'one_tail': 0.906}, 0.25: {'two_tail': 1.273, 'one_tail': 1.007}, 
            0.20: {'two_tail': 1.440, 'one_tail': 1.134}, 0.15: {'two_tail': 1.688, 'one_tail': 1.340},
            0.10: {'two_tail': 1.943, 'one_tail': 1.440}, 0.05: {'two_tail': 2.447, 'one_tail': 1.943}, 
            0.01: {'two_tail': 3.707, 'one_tail': 3.143}
        },
        7: {
            0.30: {'two_tail': 1.119, 'one_tail': 0.894}, 0.25: {'two_tail': 1.254, 'one_tail': 0.995}, 
            0.20: {'two_tail': 1.415, 'one_tail': 1.119}, 0.15: {'two_tail': 1.656, 'one_tail': 1.319},
            0.10: {'two_tail': 1.895, 'one_tail': 1.415}, 0.05: {'two_tail': 2.365, 'one_tail': 1.895}, 
            0.01: {'two_tail': 3.499, 'one_tail': 2.998}
        },
        8: {
            0.30: {'two_tail': 1.108, 'one_tail': 0.887}, 0.25: {'two_tail': 1.240, 'one_tail': 0.988}, 
            0.20: {'two_tail': 1.397, 'one_tail': 1.108}, 0.15: {'two_tail': 1.631, 'one_tail': 1.306},
            0.10: {'two_tail': 1.860, 'one_tail': 1.397}, 0.05: {'two_tail': 2.306, 'one_tail': 1.860}, 
            0.01: {'two_tail': 3.355, 'one_tail': 2.896}
        },
        9: {
            0.30: {'two_tail': 1.100, 'one_tail': 0.883}, 0.25: {'two_tail': 1.229, 'one_tail': 0.982}, 
            0.20: {'two_tail': 1.383, 'one_tail': 1.100}, 0.15: {'two_tail': 1.613, 'one_tail': 1.294},
            0.10: {'two_tail': 1.833, 'one_tail': 1.383}, 0.05: {'two_tail': 2.262, 'one_tail': 1.833}, 
            0.01: {'two_tail': 3.250, 'one_tail': 2.821}
        },
        10: {
            0.30: {'two_tail': 1.093, 'one_tail': 0.879}, 0.25: {'two_tail': 1.218, 'one_tail': 0.978}, 
            0.20: {'two_tail': 1.372, 'one_tail': 1.093}, 0.15: {'two_tail': 1.598, 'one_tail': 1.284},
            0.10: {'two_tail': 1.812, 'one_tail': 1.372}, 0.05: {'two_tail': 2.228, 'one_tail': 1.812}, 
            0.01: {'two_tail': 3.169, 'one_tail': 2.764}
        },
        12: {
            0.30: {'two_tail': 1.083, 'one_tail': 0.873}, 0.25: {'two_tail': 1.201, 'one_tail': 0.969}, 
            0.20: {'two_tail': 1.356, 'one_tail': 1.083}, 0.15: {'two_tail': 1.571, 'one_tail': 1.271},
            0.10: {'two_tail': 1.782, 'one_tail': 1.356}, 0.05: {'two_tail': 2.179, 'one_tail': 1.782}, 
            0.01: {'two_tail': 3.055, 'one_tail': 2.681}
        },
        15: {
            0.30: {'two_tail': 1.074, 'one_tail': 0.866}, 0.25: {'two_tail': 1.189, 'one_tail': 0.961}, 
            0.20: {'two_tail': 1.341, 'one_tail': 1.074}, 0.15: {'two_tail': 1.554, 'one_tail': 1.261},
            0.10: {'two_tail': 1.753, 'one_tail': 1.341}, 0.05: {'two_tail': 2.131, 'one_tail': 1.753}, 
            0.01: {'two_tail': 2.947, 'one_tail': 2.602}
        },
        20: {
            0.30: {'two_tail': 1.064, 'one_tail': 0.860}, 0.25: {'two_tail': 1.176, 'one_tail': 0.954}, 
            0.20: {'two_tail': 1.325, 'one_tail': 1.064}, 0.15: {'two_tail': 1.537, 'one_tail': 1.250},
            0.10: {'two_tail': 1.725, 'one_tail': 1.325}, 0.05: {'two_tail': 2.086, 'one_tail': 1.725}, 
            0.01: {'two_tail': 2.845, 'one_tail': 2.528}
        },
        30: {
            0.30: {'two_tail': 1.055, 'one_tail': 0.854}, 0.25: {'two_tail': 1.164, 'one_tail': 0.947}, 
            0.20: {'two_tail': 1.310, 'one_tail': 1.055}, 0.15: {'two_tail': 1.517, 'one_tail': 1.241},
            0.10: {'two_tail': 1.697, 'one_tail': 1.310}, 0.05: {'two_tail': 2.042, 'one_tail': 1.697}, 
            0.01: {'two_tail': 2.750, 'one_tail': 2.457}
        },
        60: {
            0.30: {'two_tail': 1.047, 'one_tail': 0.848}, 0.25: {'two_tail': 1.150, 'one_tail': 0.938}, 
            0.20: {'two_tail': 1.296, 'one_tail': 1.047}, 0.15: {'two_tail': 1.494, 'one_tail': 1.230},
            0.10: {'two_tail': 1.671, 'one_tail': 1.296}, 0.05: {'two_tail': 2.000, 'one_tail': 1.671}, 
            0.01: {'two_tail': 2.660, 'one_tail': 2.390}
        },
        100: {
            0.30: {'two_tail': 1.042, 'one_tail': 0.845}, 0.25: {'two_tail': 1.144, 'one_tail': 0.933}, 
            0.20: {'two_tail': 1.290, 'one_tail': 1.042}, 0.15: {'two_tail': 1.484, 'one_tail': 1.225},
            0.10: {'two_tail': 1.660, 'one_tail': 1.290}, 0.05: {'two_tail': 1.984, 'one_tail': 1.660}, 
            0.01: {'two_tail': 2.626, 'one_tail': 2.364}
        },
        'max': { # Approximates Z-distribution
            0.30: {'two_tail': 1.036, 'one_tail': 0.524}, 0.25: {'two_tail': 1.150, 'one_tail': 0.674},
            0.20: {'two_tail': 1.282, 'one_tail': 0.842}, 0.15: {'two_tail': 1.440, 'one_tail': 1.036},
            0.10: {'two_tail': 1.645, 'one_tail': 1.282}, 0.05: {'two_tail': 1.960, 'one_tail': 1.645},
            0.01: {'two_tail': 2.576, 'one_tail': 2.326}
        }
    }
    
    _CHI2_CRITICAL_VALUES = {
    # Degrees of Freedom (df):
    1: {0.30: 1.074, 0.25: 1.323, 0.20: 1.642, 0.15: 2.190, 0.10: 2.706, 0.05: 3.841, 0.01: 6.635},
    2: {0.30: 2.408, 0.25: 2.773, 0.20: 3.219, 0.15: 3.908, 0.10: 4.605, 0.05: 5.991, 0.01: 9.210},
    3: {0.30: 3.665, 0.25: 4.108, 0.20: 4.642, 0.15: 5.518, 0.10: 6.251, 0.05: 7.815, 0.01: 11.345},
    4: {0.30: 4.878, 0.25: 5.437, 0.20: 6.064, 0.15: 7.070, 0.10: 7.779, 0.05: 9.488, 0.01: 13.277},
    5: {0.30: 6.064, 0.25: 6.626, 0.20: 7.289, 0.15: 8.328, 0.10: 9.236, 0.05: 11.070, 0.01: 15.086},
    6: {0.30: 7.231, 0.25: 7.964, 0.20: 8.558, 0.15: 9.774, 0.10: 10.645, 0.05: 12.592, 0.01: 16.812},
    7: {0.30: 8.384, 0.25: 9.037, 0.20: 9.803, 0.15: 11.205, 0.10: 12.017, 0.05: 14.067, 0.01: 18.475},
    8: {0.30: 9.524, 0.25: 10.151, 0.20: 11.030, 0.15: 12.442, 0.10: 13.361, 0.05: 15.507, 0.01: 20.090},
    9: {0.30: 10.656, 0.25: 11.261, 0.20: 12.242, 0.15: 13.719, 0.10: 14.684, 0.05: 16.919, 0.01: 21.666},
    10: {0.30: 11.781, 0.25: 12.382, 0.20: 13.442, 0.15: 14.980, 0.10: 15.987, 0.05: 18.307, 0.01: 23.209},
    12: {0.30: 14.001, 0.25: 14.707, 0.20: 16.000, 0.15: 17.581, 0.10: 18.549, 0.05: 21.026, 0.01: 26.217},
    15: {0.30: 17.653, 0.25: 18.477, 0.20: 19.960, 0.15: 22.307, 0.10: 22.307, 0.05: 24.996, 0.01: 30.578},
    20: {0.30: 22.956, 0.25: 24.433, 0.20: 25.038, 0.15: 28.412, 0.10: 28.412, 0.05: 31.410, 0.01: 37.566},
    30: {0.30: 34.800, 0.25: 36.452, 0.20: 38.300, 0.15: 40.256, 0.10: 40.256, 0.05: 43.773, 0.01: 50.892},
    60: {0.30: 65.410, 0.25: 68.300, 0.20: 70.000, 0.15: 72.825, 0.10: 74.397, 0.05: 79.082, 0.01: 88.379},
    100: {0.30: 108.611, 0.25: 112.500, 0.20: 115.000, 0.15: 118.498, 0.10: 118.498, 0.05: 124.342, 0.01: 135.807},
    # Chi-square tables typically don't use a 'max' Z-approximation 
    # as the distribution shape is different, but you can approximate
    # using the formula for large df: chi2 = df + z * sqrt(2*df)
    }   
    
    def __repr__(self):
        """If the toolbox is printed, display a message."""
        return ".stats: a set of functions for performing statistical operations on lists of numerical values."

    def _check_significance(self, stat_value, df, test_type, alpha=0.05, tail_type='two_tail'):
        """
        Retrieves the critical value from the internal tables and compares it
        to the calculated test statistic to determine statistical significance.

        Degrees of freedom (df) is rounded down to the nearest available table value
        (e.g., df=13 rounds to 12) to ensure a conservative test.

        :param stat_value: The calculated T or Chi-Square statistic.
        :param df: Degrees of Freedom (integer).
        :param test_type: 't' for t-test or 'chi2' for chi-square.
        :param alpha: The significance level (float, e.g., 0.05).
        :param tail_type: For t-tests, 'one_tail' or 'two_tail'. Ignored for chi-square.
        :return: A string indicating significance and the critical value.
        """
        test_type = test_type.lower()
        
        # 1. Select the correct dictionary
        if test_type == 't':
            crit_vals_dict = self._T_CRITICAL_VALUES
        elif test_type == 'chi2':
            crit_vals_dict = self._CHI2_CRITICAL_VALUES
        else:
            return f"Error: Unknown test_type '{test_type}'. Must be 't' or 'chi2'."

        # 2. Get and sort the integer DF keys for lookup
        df_keys = [key for key in crit_vals_dict.keys() if isinstance(key, int)]
        df_keys.sort() # e.g., [1, 2, 3, ..., 60, 100]
        
        # 3. Determine the closest *lower* degrees of freedom key (Conservative Rounding)
        df_key = 1 # Default to 1 if df is extremely small
        
        if df >= df_keys[-1]:
            # Use 'max' for T-test if >= 100, otherwise use 100 for Chi2
            df_key = 'max' if (test_type == 't' and 'max' in crit_vals_dict) else df_keys[-1]
        elif df < df_keys[0]:
            df_key = df_keys[0] # Use lowest DF (1)
        else:
            # Find the largest key smaller than df
            for key in df_keys:
                if key <= df:
                    df_key = key
                else:
                    break # Since the list is sorted, the last one found is the correct one

        # 4. Retrieve the critical value and prepare the statistic
        try:
            if test_type == 't':
                # T-tests use the absolute value of the statistic
                stat_for_comp = abs(stat_value)
                critical_val = crit_vals_dict[df_key][alpha][tail_type]
            else: # chi2
                # Chi-Square tests are one-tailed and inherently positive
                stat_for_comp = stat_value
                critical_val = crit_vals_dict[df_key][alpha]
        except KeyError:
            return f"Error: Failed to find critical value for df={df_key} and alpha={alpha}. Check dictionary structure."

        # 5. Compare and return result
        # Significance is achieved when the calculated statistic exceeds the critical value.
        if stat_for_comp > critical_val:
            result_string = "Likely Significant (Reject H0)"
        else:
            result_string = "Likely Not Significant (Fail to Reject H0)"
        
        return f"{result_string}. Critical Value: {critical_val} (DF used: {df_key})."

    def display_toolbox_functions(self):
        """Display a list of all available functions within this toolbox."""
        print(f"Number of {__class__.__name__[1:]} functions: {len([f for f in dir(__class__) if not f.startswith('_')])}")
        for f in [f for f in dir(__class__) if not f.startswith("_")]:
            print(f)

    # ============================================================
    # DESCRIPTIVE STATISTICS
    # ============================================================

    def descriptive_stat(self,vals,stat,decimals = 2):
        """Calculate and return the specified descriptive statistic ('stat') of the input list ('vals').

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        stat (string): the specific statistic to return. Possible statistics are "Count", "Unique Values", "First", "Last", "Sum", "Max", "Min", "Range", "Mean", "Median", "Mode", "IQ1", "IQ3", "IQR", "Standard Deviation", "Standard Error", "Variance", "Coefficient of Variation", "Skewness", and "Kurtosis"
        
        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        valid_stats = {
            "Count", "Unique Values", "First", "Last", "Sum", "Max", "Min", "Range", "Mean", "Median", "Mode",
            "IQ1", "IQ3", "IQR", "Standard Deviation", "Standard Error", "Variance",
            "Coefficient of Variation", "Skewness", "Kurtosis"
        }
        if not isinstance(stat, str):
            raise TypeError("'stat' must be a string.")
        if stat not in valid_stats:
            raise ValueError(f"Statistic '{stat}' is not recognized.")
        # Basic stats
        if stat == "Count":
            return len(vals)
        if stat == "Unique Values":
            return len(set(vals))
        if stat == "First":
            return vals[0]
        if stat == "Last":
            return vals[-1]
        if stat == "Sum":
            return round(sum(vals), decimals)
        if stat == "Max":
            return round(max(vals), decimals)
        if stat == "Min":
            return round(min(vals), decimals)
        if stat == "Range":
            return round(max(vals) - min(vals), decimals)
        if stat == "Mean":
            return round(sum(vals) / len(vals), decimals)
        if stat == "Median":
            return round(statistics.median(vals), decimals)
        if stat == "Mode":
            try:
                return round(statistics.mode(vals), decimals)
            except:
                return "None"
        # Quartiles and IQR
        lv = len(vals)
        s_v = sorted(vals)
        if stat == "IQ1":
            if lv % 2 == 0:
                return round(statistics.median(s_v[:lv // 2]), decimals)
            else:
                mid = (lv - 1) // 2
                return round(statistics.median(s_v[:mid]), decimals)
        if stat == "IQ3":
            if lv % 2 == 0:
                return round(statistics.median(s_v[lv // 2:]), decimals)
            else:
                mid = (lv - 1) // 2
                return round(statistics.median(s_v[mid + 1:]), decimals)
        if stat == "IQR":
            if lv % 2 == 0:
                q1 = statistics.median(s_v[:lv // 2])
                q3 = statistics.median(s_v[lv // 2:])
            else:
                mid = (lv - 1) // 2
                q1 = statistics.median(s_v[:mid])
                q3 = statistics.median(s_v[mid + 1:])
            return round(q3 - q1, decimals)
        # Higher-order stats
        mean = sum(vals) / lv
        sum_sqr = 0
        sum_cube = 0
        sum_quart = 0
        for v in vals:
            diff = v - mean
            sum_sqr += diff ** 2
            sum_cube += diff ** 3
            sum_quart += diff ** 4
        if stat == "Standard Deviation":
            return round(math.sqrt(sum_sqr / lv), decimals)
        if stat == "Standard Error":
            std_dev = math.sqrt(sum_sqr / lv)
            return round(std_dev / math.sqrt(lv), decimals)
        if stat == "Variance":
            std_dev = math.sqrt(sum_sqr / lv)
            return round(std_dev * std_dev, decimals)
        if stat == "Coefficient of Variation":
            std_dev = math.sqrt(sum_sqr / lv)
            return round(std_dev / mean, decimals)
        if stat == "Skewness":
            variance = sum_sqr / lv
            return round((sum_cube / lv) / (variance ** 1.5), decimals)
        if stat == "Kurtosis":
            std_dev = math.sqrt(sum_sqr / lv)
            return round((sum_quart / lv) / (std_dev ** 4), decimals)

    def summary_stats(self, vals, terminal = False, decimals=2):
        """Calculate and display all descriptive statistics of the input list ('vals') and return a dictionary of all statistics.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        terminal (Boolean): flag to indicate if the result will be printed to a terminal or just returned.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        from textwrap import dedent
        lv = len(vals)
        vsum = round(sum(vals), decimals)
        mean = round(vsum / lv, decimals)
        vmax = round(max(vals), decimals)
        vmin = round(min(vals), decimals)
        data_range = round(vmax - vmin, decimals)
        medi = round(statistics.median(vals), decimals)
        try:
            mode = round(statistics.mode(vals), decimals)
        except:
            mode = "None"
        # Higher-order sums
        sum_sqr = 0
        sum_cube = 0
        sum_quart = 0
        unique_vals = []
        for v in vals:
            if v not in unique_vals:
                unique_vals.append(v)
            diff = v - mean
            sum_sqr += diff ** 2
            sum_cube += diff ** 3
            sum_quart += diff ** 4
        # Quartiles and IQR
        s_v = sorted(vals)
        if lv % 2 == 0:
            q1 = statistics.median(s_v[:lv // 2])
            q3 = statistics.median(s_v[lv // 2:])
        else:
            mid = (lv - 1) // 2
            q1 = statistics.median(s_v[:mid])
            q3 = statistics.median(s_v[mid + 1:])
        q1 = round(q1, decimals)
        q3 = round(q3, decimals)
        iqr = round(q3 - q1, decimals)
        # Standard deviation, variance, etc.
        variance = sum_sqr / lv
        std_dev = round(math.sqrt(variance), decimals)
        std_err = round(std_dev / math.sqrt(lv), decimals)
        var = round(variance, decimals)
        cvar = round(std_dev / mean, decimals)
        # Skewness and kurtosis (corrected formulas)
        skew = round((sum_cube / lv) / (variance ** 1.5), decimals)
        kurt = round((sum_quart / lv) / (variance ** 2), decimals)
        output = f"""
        Count: {lv}
        Unique Values: {len(unique_vals)}
        First: {vals[0]}
        Last: {vals[-1]}
        Sum: {vsum}
        Max: {vmax}
        Min: {vmin}
        Range: {data_range}
        Mean: {mean}
        Median: {medi}
        Mode: {mode}
        Lower Quartile (IQ1): {q1}
        Upper Quartile (IQ3): {q3}
        Interquartile Range (IQR): {iqr}
        Standard Deviation: {std_dev}
        Standard Error: {std_err}
        Variance: {var}
        Coefficient of Variation: {cvar}
        Skewness: {skew}
        Kurtosis: {kurt}
        """
        if terminal:
            print(dedent(output).lstrip())
        return {
            "Count": lv,
            "Unique Values": len(unique_vals),
            "First": vals[0],
            "Last": vals[-1],
            "Sum": vsum,
            "Max": vmax,
            "Min": vmin,
            "Range": data_range,
            "Mean": mean,
            "Median": medi,
            "Mode": mode,
            "IQ1": q1,
            "IQ3": q3,
            "IQR": iqr,
            "Standard Deviation": std_dev,
            "Standard Error": std_err,
            "Variance": var,
            "Coefficient of Variation": cvar,
            "Skewness": skew,
            "Kurtosis": kurt
        }

    def root_mean_square_error(self,predicted,observed,decimals = 2):
        """Calculate and return the root-mean-square-error for input lists of predicted ('predicted') and observed ('observed') values.
        
        Parameters:
        
        predicted (list): the list of predicted values the operation will be performed on. The input list will not be altered by the operation.

        observed (list): the list of observed values the operation will be performed on. The input list will not be altered by the operation.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(predicted, list) or not predicted:
            raise TypeError("'predicted' must be a non-empty list.")
        if not isinstance(observed, list) or not observed:
            raise TypeError("'observed' must be a non-empty list.")

        if len(predicted) != len(observed):
            raise ValueError("'predicted' and 'observed' must be the same length.")

        if not all(isinstance(v, (int, float)) for v in predicted):
            raise ValueError("All elements in 'predicted' must be numeric.")
        if not all(isinstance(v, (int, float)) for v in observed):
            raise ValueError("All elements in 'observed' must be numeric.")

        n = len(observed)
        rsum = 0

        for i in range(n):
            rsum += (predicted[i] - observed[i]) ** 2

        return round(math.sqrt(rsum / n), decimals)
    
    def mean_square_error(self,predicted,observed,decimals = 2):
        """Calculate and return the mean square error for input lists of predicted ('predicted') and observed ('observed') values.
        
        Parameters:
        
        predicted (list): the list of predicted values the operation will be performed on. The input list will not be altered by the operation.

        observed (list): the list of observed values the operation will be performed on. The input list will not be altered by the operation.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(predicted, list) or not predicted:
            raise TypeError("'predicted' must be a non-empty list.")
        if not isinstance(observed, list) or not observed:
            raise TypeError("'observed' must be a non-empty list.")

        if len(predicted) != len(observed):
            raise ValueError("'predicted' and 'observed' must be the same length.")

        if not all(isinstance(v, (int, float)) for v in predicted):
            raise ValueError("All elements in 'predicted' must be numeric.")
        if not all(isinstance(v, (int, float)) for v in observed):
            raise ValueError("All elements in 'observed' must be numeric.")

        n = len(observed)
        rsum = 0

        for i in range(n):
            rsum += (predicted[i] - observed[i]) ** 2

        return round(rsum / n, decimals)

    def mean_absolute_error(self,predicted,observed,decimals = 2):
        """Calculate and return the mean absolute error for input lists of predicted ('predicted') and observed ('observed') values.
        
        Parameters:
        
        predicted (list): the list of predicted values the operation will be performed on. The input list will not be altered by the operation.

        observed (list): the list of observed values the operation will be performed on. The input list will not be altered by the operation.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(predicted, list) or not predicted:
            raise TypeError("'predicted' must be a non-empty list.")
        if not isinstance(observed, list) or not observed:
            raise TypeError("'observed' must be a non-empty list.")

        if len(predicted) != len(observed):
            raise ValueError("'predicted' and 'observed' must be the same length.")

        if not all(isinstance(v, (int, float)) for v in predicted):
            raise ValueError("All elements in 'predicted' must be numeric.")
        if not all(isinstance(v, (int, float)) for v in observed):
            raise ValueError("All elements in 'observed' must be numeric.")

        n = len(observed)
        rsum = 0

        for i in range(n):
            rsum += abs(predicted[i] - observed[i])

        return round(rsum / n, decimals)

    def single_freq(self,vals,value):
        """Calculate and return the frequency with which a specified value ('value') appears in the input list ('vals').

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        value (integer, float, string): the value to determine the frequency of."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        return vals.count(value)
       
    def unique_freq(self,vals):
        """Determine the unique values in the input list ('vals') and calculate the frequency with which each unique value appears. Two lists are returned: a list of unique values and a list of unique value frequencies.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        counts = {}
        order = []

        for v in vals:
            if v in counts:
                counts[v] += 1
            else:
                counts[v] = 1
                order.append(v)

        # Build parallel lists in order of first appearance
        unique_vals = order
        freq = [counts[v] for v in order]

        return unique_vals, freq

    def k_most_common(self, vals, k):
        """Determine the unique values in the input list ('vals') and return the k most frequent values.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        k (integer): the number of most frequent values to return."""
        unique_vals, counts = self.unique_freq(vals)
        freq = sorted(zip(unique_vals, counts), key=lambda x: x[1], reverse=True)
        return [val for val, _ in freq[:k]]

    # ============================================================
    # CONTINGENCY ANALYSIS
    # ============================================================

    def contingency_table(self, x_var, y_var, sort=True):
        """Construct a contingency table (cross-tabulation) from two categorical variables of equal length. Returns a dictionary containing the contingency table (key: c_table), row labels (key: row_labels), and column labels (key: col_labels).

        Parameters:

        x_var (list): the list of categorical values representing the first variable. The input list will not be altered by the operation.

        y_var (list): the list of categorical values representing the second variable. The input list will not be altered by the operation.

        sort (bool): if True, categories are sorted alphabetically. If False, they appear in the order of first occurrence."""
        # Validation
        if not isinstance(x_var, list) or not x_var:
            raise TypeError("'x_var' must be a non-empty list.")
        if not isinstance(y_var, list) or not y_var:
            raise TypeError("'y_var' must be a non-empty list.")
        if len(x_var) != len(y_var):
            raise ValueError("'x_var' and 'y_var' must be the same length.")
        if not all(isinstance(v, (str, int, float)) for v in x_var):
            raise ValueError("All elements in 'x_var' must be categorical (string or number).")
        if not all(isinstance(v, (str, int, float)) for v in y_var):
            raise ValueError("All elements in 'y_var' must be categorical (string or number).")

        # Identify unique categories
        if sort:
            row_labels = sorted(set(x_var))
            col_labels = sorted(set(y_var))
        else:
            # Preserve first occurrence order
            row_labels = []
            col_labels = []
            for v in x_var:
                if v not in row_labels:
                    row_labels.append(v)
            for v in y_var:
                if v not in col_labels:
                    col_labels.append(v)

        # Initialize contingency matrix
        c_table = [[0 for _ in col_labels] for _ in row_labels]

        # Build index maps for fast lookup
        row_index = {label: i for i, label in enumerate(row_labels)}
        col_index = {label: i for i, label in enumerate(col_labels)}

        # Count co-occurrences
        for xv, yv in zip(x_var, y_var):
            r = row_index[xv]
            c = col_index[yv]
            c_table[r][c] += 1

        return {
            "c_table": c_table,
            "row_labels": row_labels,
            "col_labels": col_labels
        }

    def phi_coefficient(self, c_table, decimals=2):
        """Calculate the Phi coefficient (effect size) for a 2x2 contingency table and return a dictionary containing phi (key: phi) and total sample size (key: n)

        Parameters:

        c_table (list): a nested list representing the numeric contingency matrix. This should be the 'c_table' value extracted from the dictionary returned by the contingency_table() function, not the full dictionary itself.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(c_table, list) or len(c_table) != 2:
            raise ValueError("Phi coefficient requires a 2x2 contingency table.")
        if not all(isinstance(row, list) and len(row) == 2 for row in c_table):
            raise ValueError("Phi coefficient requires a 2x2 contingency table.")
        if not all(all(isinstance(v, (int, float)) and v >= 0 for v in row) for row in c_table):
            raise ValueError("All values in the contingency table must be non-negative numbers.")

        # Total sample size
        n = sum(sum(row) for row in c_table)

        # Compute phi using the determinant method
        a, b = c_table[0]
        c, d = c_table[1]

        numerator = (a*d - b*c) ** 2 * n
        denominator = (a + b) * (c + d) * (a + c) * (b + d)

        if denominator == 0:
            raise ValueError("Invalid contingency table: division by zero in phi calculation.")

        phi = math.sqrt(numerator / denominator)
        phi = round(phi, decimals)

        return {
            "phi": phi,
            "n": n
        }
    
    def cramers_v(self, c_table, decimals=2):
        """Calculate Cramér's V (effect size) for any r x c contingency table and return a dictionary containing V (key: cramers_v), chi-square (key: chi_square), degrees of freedom (key: df), and sample size (key: n).

        Parameters:

        c_table (list): a nested list representing the numeric contingency matrix. This should be the 'c_table' value extracted from the dictionary returned by the contingency_table() function, not the full dictionary itself.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(c_table, list) or not c_table:
            raise TypeError("'c_table' must be a non-empty list of lists.")
        if not all(isinstance(row, list) and row for row in c_table):
            raise TypeError("Each row in 'c_table' must be a non-empty list.")
        if len({len(row) for row in c_table}) != 1:
            raise ValueError("All rows in the contingency table must have the same length.")
        if not all(all(isinstance(v, (int, float)) and v >= 0 for v in row) for row in c_table):
            raise ValueError("All values in the contingency table must be non-negative numbers.")

        # Use your chi-square independence function
        chi_result = self.chi_square_independence(c_table, alpha=0.05, decimals=decimals)

        chi_sqr = chi_result["chi_square"]
        df = chi_result["df"]
        n = sum(sum(row) for row in c_table)

        r = len(c_table)
        c = len(c_table[0])
        k = min(r, c)

        if k <= 1:
            raise ValueError("Cramér's V requires at least a 2x2 table.")

        V = math.sqrt(chi_sqr / (n * (k - 1)))
        V = round(V, decimals)

        return {
            "cramers_v": V,
            "chi_square": chi_sqr,
            "df": df,
            "n": n
        }

    def binary_stats(self, predicted, observed, beta=1, decimals=2):
        """Calculate and return binary classification statistics for lists of predicted ('predicted') and observed ('observed') values. Returns a dictionary with the following keys: recall, precision, accuracy, f1_score, phi_coefficient, confusion_matrix (a dictionary of TP, FP, TN, FN), fbeta_score, balanced_accuracy, negative_predictive_value, specificity, and kappa.

        Parameters:

        predicted (list): the list of binary values (1/0 or True/False) representing the data produced by a classification technique. The input list will not be altered by the operation.

        observed (list): the list of binary values (1/0 or True/False) representing the observed data to compare the classification technique against. The input list will not be altered by the operation.

        beta (float): the beta parameter for the F-beta score (default is 1 for F1).

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(predicted, list) or not predicted:
            raise TypeError("'predicted' must be a non-empty list.")
        if not isinstance(observed, list) or not observed:
            raise TypeError("'observed' must be a non-empty list.")
        if len(predicted) != len(observed):
            raise ValueError("'predicted' and 'observed' must be the same length.")

        # Ensure binary values
        valid = {0, 1, True, False}
        if not all(v in valid for v in predicted):
            raise ValueError("'predicted' must contain only binary values (1/0 or True/False).")
        if not all(v in valid for v in observed):
            raise ValueError("'observed' must contain only binary values (1/0 or True/False).")

        tp = tn = fp = fn = 0

        for i in range(len(predicted)):
            p = bool(predicted[i])
            o = bool(observed[i])

            if p and o:
                tp += 1
            elif p and not o:
                fp += 1
            elif not p and o:
                fn += 1
            else:
                tn += 1

        # Avoid division by zero
        recall = tp / (tp + fn) if (tp + fn) else 0
        precision = tp / (tp + fp) if (tp + fp) else 0
        accuracy = (tp + tn) / (tp + tn + fp + fn)

        # F1 score
        f_score = (2 * tp) / ((2 * tp) + fp + fn) if ((2 * tp) + fp + fn) else 0

        # F-beta score
        beta_sq = beta * beta
        fbeta_num = (1 + beta_sq) * tp
        fbeta_den = (1 + beta_sq) * tp + beta_sq * fn + fp
        fbeta_score = fbeta_num / fbeta_den if fbeta_den else 0

        # Matthews Correlation Coefficient (Pearson's Phi)
        denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        phi = ((tp * tn) - (fp * fn)) / denom if denom else 0

        # Specificity (True Negative Rate)
        specificity = tn / (tn + fp) if (tn + fp) else 0

        # Negative Predictive Value
        npv = tn / (tn + fn) if (tn + fn) else 0

        # Balanced Accuracy
        balanced_accuracy = (recall + specificity) / 2

        # Cohen's Kappa
        total = tp + tn + fp + fn
        po = accuracy  # observed agreement
        pe = (((tp + fp) * (tp + fn)) + ((fn + tn) * (fp + tn))) / (total * total)
        kappa = (po - pe) / (1 - pe) if (1 - pe) else 0

        # Confusion Matrix dictionary
        confusion_matrix = {
            "TP": tp,
            "FP": fp,
            "TN": tn,
            "FN": fn
        }

        return {
            "recall": round(recall, decimals),
            "precision": round(precision, decimals),
            "accuracy": round(accuracy, decimals),
            "f1_score": round(f_score, decimals),
            "phi_coefficient": round(phi, decimals),
            "fbeta_score": round(fbeta_score, decimals),
            "specificity": round(specificity, decimals),
            "negative_predictive_value": round(npv, decimals),
            "balanced_accuracy": round(balanced_accuracy, decimals),
            "kappa": round(kappa, decimals),
            "confusion_matrix": confusion_matrix
        }

    # ============================================================
    # DISTRIBUTION AND DENSITY FUNCTIONS
    # ============================================================
    
    def probability_distribution(self, vals):
        """Calculate the probability distribution for an input list (vals). Returns a dictionary containing the unique values (key: values), their counts (key: counts), their probabilities (key: probabilities), and the total number of observations (key: total).

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation."""
        # Validation
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")

        # Count occurrences
        counts = {}
        for v in vals:
            counts[v] = counts.get(v, 0) + 1

        total = len(vals)

        values = list(counts.keys())
        count_list = [counts[v] for v in values]
        probabilities = [c / total for c in count_list]

        return {
            "values": values,
            "counts": count_list,
            "probabilities": probabilities,
            "total": total
        }

    def cumulative_mass_function(self, vals, decimals=2):
        """Calculate the cumulative mass function for an input list ('vals)'. Returns a dictionary containing the cumulative values (key: cumulative_values), the normalized cumulative mass function (key: cmf), and the total sum of values (key: total).

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")

        total = sum(vals)
        if total == 0:
            raise ValueError("The sum of 'vals' must not be zero.")

        cumulative_values = []
        cmf = []

        running_total = 0
        for v in vals:
            running_total += v
            cumulative_values.append(running_total)
            cmf.append(round(running_total / total, decimals))

        return {
            "cumulative_values": cumulative_values,
            "cmf": cmf,
            "total": total
        }

    def cumulative_distribution(self, vals, decimals=2):
        """Calculate the cumulative distribution function for an input list ('vals'). Returns a dictionary containing the sorted values (key: values), their probabilities (key: probabilities), the cumulative distribution (key: cdf), and the total number of observations (key: total).

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")

        total = len(vals)

        # Sort values
        sorted_vals = sorted(vals)

        # Compute probabilities (each value has probability 1/n)
        prob = 1 / total
        probabilities = [round(prob, decimals) for _ in sorted_vals]

        # Compute cumulative distribution
        cdf = []
        running = 0
        for _ in sorted_vals:
            running += prob
            cdf.append(round(running, decimals))

        return {
            "values": sorted_vals,
            "probabilities": probabilities,
            "cdf": cdf,
            "total": total
        }

    # ============================================================
    # STANDARDIZATION, SAMPLING, CLEANING
    # ============================================================

    def get_random_sample(self,vals,num_samples):
        """Return a new list containing a specified number of random samples ('num_samples') from the input list ('vals'), without replacement.

        Parameters:

        vals (list): the list of values to draw random samples from. The input list will not be altered by the operation.

        num_samples (integer): the number of random samples to draw from the given list."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not isinstance(num_samples, int) or num_samples <= 0:
            raise ValueError("'num_samples' must be a positive integer.")
        if num_samples > len(vals):
            raise ValueError("'num_samples' cannot be greater than the length of 'vals'.")
        return random.sample(vals,num_samples)
    
    def remove_beyond_std(self,vals,std_devs,decimals = 2):
        """Remove values from an input list ('vals') that fall outside of the specified number of standard deviations (std_devs) and return a new list. Can be used for outlier removal.


        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        std_dev (integer): the number of standard deviations to remove values by.
        
        decimals (integer): the number of decimal places the mean and standard deviation will be calculated to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")

        if not isinstance(std_devs, (int, float)) or std_devs <= 0:
            raise ValueError("'std_devs' must be a positive integer or float.")

        # Compute mean and standard deviation using your validated descriptive_stat
        mean = self.descriptive_stat(vals, "Mean", decimals)
        std = self.descriptive_stat(vals, "Standard Deviation", decimals)

        thresh_upper = mean + (std * std_devs)
        thresh_lower = mean - (std * std_devs)

        output = []
        for v in vals:
            if v > thresh_lower and v < thresh_upper:
                output.append(v)

        return output

    def remove_beyond_iqr(self,vals):
        """Remove values from an input list ('vals') that fall outside of the interquartile range and return a new list. Can be used for outlier removal.
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")

        lv = len(vals)
        s_v = sorted(vals)

        # Compute Q1 and Q3 using the same logic as descriptive_stat
        if lv % 2 == 0:
            q1 = statistics.median(s_v[:lv // 2])
            q3 = statistics.median(s_v[lv // 2:])
        else:
            mid = (lv - 1) // 2
            q1 = statistics.median(s_v[:mid])
            q3 = statistics.median(s_v[mid + 1:])

        iqr = q3 - q1

        lower = q1 - 1.5 * iqr
        upper = q3 + 1.5 * iqr

        output = []
        for v in vals:
            if v > lower and v < upper:
                output.append(v)

        return output

    def standardize(self,vals,decimals = 2):
        """Calculate the z-score, or standard score, of each value in the input list ('vals') and return a new list.
        
        Parameters:
        
        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")

        lv = len(vals)
        mean = sum(vals) / lv

        # Compute variance (population)
        sum_sqr = 0
        for v in vals:
            sum_sqr += (v - mean) ** 2

        variance = sum_sqr / lv
        std_dev = math.sqrt(variance)

        if std_dev == 0:
            raise ValueError("Standard deviation is zero; z-scores cannot be computed.")

        return [round((v - mean) / std_dev, decimals) for v in vals]
    
    # ============================================================
    # CORRELATION
    # ============================================================

    def pearson_correlation(self, x_var, y_var, decimals=2):
        """Calculate and return Pearson's correlation coefficient between an x variable (x_var) and a y variable (y_var). Returns a dictionary containing the correlation coefficient (key: correlation) and covariance (key: covariance).

        Parameters:

        x_var (list): the list of values representing the x variable the operation will be performed on. The input list will not be altered by the operation.

        y_var (list): the list of values representing the y variable the operation will be performed on. The input list will not be altered by the operation.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(x_var, list) or not x_var:
            raise TypeError("'x_var' must be a non-empty list.")
        if not isinstance(y_var, list) or not y_var:
            raise TypeError("'y_var' must be a non-empty list.")
        if len(x_var) != len(y_var):
            raise ValueError("'x_var' and 'y_var' must be the same length.")
        if not all(isinstance(v, (int, float)) for v in x_var):
            raise ValueError("All elements in 'x_var' must be numeric.")
        if not all(isinstance(v, (int, float)) for v in y_var):
            raise ValueError("All elements in 'y_var' must be numeric.")

        n = len(x_var)
        xmean = sum(x_var) / n
        ymean = sum(y_var) / n

        covar = 0
        sum_sqr_x = 0
        sum_sqr_y = 0

        for i in range(n):
            dx = x_var[i] - xmean
            dy = y_var[i] - ymean
            covar += dx * dy
            sum_sqr_x += dx ** 2
            sum_sqr_y += dy ** 2

        # Pearson correlation coefficient
        correlation = covar / (math.sqrt(sum_sqr_x) * math.sqrt(sum_sqr_y))

        # Sample covariance
        covariance = covar / (n - 1)

        return {
            "correlation": round(correlation, decimals),
            "covariance": round(covariance, decimals)
        }
    
    def spearman_rank_correlation(self,x_var,y_var,decimals = 2):
        """Calculate and return the Spearman Rank-Order Correlation Coefficient between an x variable ('x_var') and a y variable ('y_var').

        Parameters:

        x_var (list): the list of values representing the ranked x variable the operation will be performed on. Values must be integers. The input list will not be altered by the operation.

        y_var (list): the list of values representing the ranked y variable the operation will be performed on. Values must be integers. The input list will not be altered by the operation.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(x_var, list) or not x_var:
            raise TypeError("'x_var' must be a non-empty list.")
        if not isinstance(y_var, list) or not y_var:
            raise TypeError("'y_var' must be a non-empty list.")
        if len(x_var) != len(y_var):
            raise ValueError("'x_var' and 'y_var' must be the same length.")

        if not all(isinstance(v, int) for v in x_var):
            raise ValueError("All elements in 'x_var' must be integers representing ranks.")
        if not all(isinstance(v, int) for v in y_var):
            raise ValueError("All elements in 'y_var' must be integers representing ranks.")

        n = len(x_var)
        rsum = 0

        for i in range(n):
            rsum += (x_var[i] - y_var[i]) ** 2

        correlation = 1 - ((6 * rsum) / (n * (n**2 - 1)))

        return round(correlation, decimals)

    def kendall_tau(self, x_var, y_var, alpha=0.05, decimals=2):
        """Compute Kendall's Tau-b correlation coefficient for two numeric variables ('x_var', 'y_var'). Returns a dictionary containing tau (key: tau), the Z statistic (key: z_stat), significance decision (key: significance), and the critical value used (key: critical_value).

        Parameters:

        x_var (list): the list of numeric values representing the first variable. The input list will not be altered by the operation.

        y_var (list): the list of numeric values representing the second variable. The input list will not be altered by the operation.

        alpha (float): the significance level for the Kendall test.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(x_var, list) or not x_var:
            raise TypeError("'x_var' must be a non-empty list.")
        if not isinstance(y_var, list) or not y_var:
            raise TypeError("'y_var' must be a non-empty list.")
        if len(x_var) != len(y_var):
            raise ValueError("'x_var' and 'y_var' must be the same length.")
        if not all(isinstance(v, (int, float)) for v in x_var):
            raise ValueError("All elements in 'x_var' must be numeric.")
        if not all(isinstance(v, (int, float)) for v in y_var):
            raise ValueError("All elements in 'y_var' must be numeric.")
        if not isinstance(alpha, float):
            raise TypeError("'alpha' must be a float value.")

        n = len(x_var)

        # Count concordant, discordant, and ties
        concordant = 0
        discordant = 0
        ties_x = 0
        ties_y = 0

        for i in range(n - 1):
            for j in range(i + 1, n):
                dx = x_var[i] - x_var[j]
                dy = y_var[i] - y_var[j]

                if dx == 0 and dy == 0:
                    # tied on both — ignored in tau-b numerator
                    continue
                elif dx == 0:
                    ties_x += 1
                elif dy == 0:
                    ties_y += 1
                else:
                    # concordant or discordant
                    if dx * dy > 0:
                        concordant += 1
                    else:
                        discordant += 1

        # Tau-b computation
        numerator = concordant - discordant
        denom_x = (n * (n - 1) / 2) - ties_x
        denom_y = (n * (n - 1) / 2) - ties_y

        if denom_x == 0 or denom_y == 0:
            raise ValueError("Kendall's tau is undefined due to excessive ties.")

        tau = numerator / (denom_x * denom_y) ** 0.5
        tau = round(tau, decimals)

        # Normal approximation for Z
        # Standard error under H0
        se = (2 * (2*n + 5)) / (9 * n * (n - 1))
        se = se ** 0.5

        z_stat = tau / se
        z_stat = round(z_stat, decimals)

        #  Significance check using Z row of t-table
        significance_result = self._check_significance(
            stat_value=z_stat,
            df=999999,  # triggers Z-distribution row
            test_type='t',
            alpha=alpha,
            tail_type='two_tail'
        )

        # Extract critical value and df_used
        try:
            crit_part = significance_result.split("Critical Value:")[1]
            critical_value = float(crit_part.split("(")[0].strip())
            df_used = crit_part.split("DF used:")[1].replace(").", "").strip()
        except Exception:
            critical_value = None
            df_used = None

        return {
            "tau": tau,
            "z_stat": z_stat,
            "n": n,
            "alpha": alpha,
            "critical_value": critical_value,
            "df_used": df_used,
            "significance": significance_result
        }

    # ============================================================
    # REGRESSION MODELLING
    # ============================================================

    def linear_regression(self,x_var,y_var,decimals = 2):
        """Calculate the equation of the regression line between an independent variable ('x_var') and a dependent variable ('y_var'). Returns a dictionary containing the intercept (key: intercept), slope (key: slope), r-squared (key: r_squared), standard error (key: std_error), and residuals (key: residuals).

        Parameters:

        x_var (list): the list of values representing the independent variable the operation will be performed on. The input list will not be altered by the operation.

        y_var (list): the list of values representing the dependent variable the operation will be performed on. The input list will not be altered by the operation.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(x_var, list) or not x_var:
            raise TypeError("'x_var' must be a non-empty list.")
        if not isinstance(y_var, list) or not y_var:
            raise TypeError("'y_var' must be a non-empty list.")
        if len(x_var) != len(y_var):
            raise ValueError("'x_var' and 'y_var' must be the same length.")
        if not all(isinstance(v, (int, float)) for v in x_var):
            raise ValueError("All elements in 'x_var' must be numeric.")
        if not all(isinstance(v, (int, float)) for v in y_var):
            raise ValueError("All elements in 'y_var' must be numeric.")
        if len(x_var) < 3:
            raise ValueError("At least 3 data points are required for regression.")

        # Means
        xmean = sum(x_var) / len(x_var)
        ymean = sum(y_var) / len(y_var)

        # Sums for slope and correlation components
        x_sqr_sum = 0
        x_sqr_sum2 = 0
        y_sqr_sum = 0

        for i in range(len(x_var)):
            x_sqr_sum += (x_var[i] - xmean) * (y_var[i] - ymean)
            x_sqr_sum2 += (x_var[i] - xmean) ** 2
            y_sqr_sum += (y_var[i] - ymean) ** 2

        # Regression coefficients
        b1 = x_sqr_sum / x_sqr_sum2
        b_nought = ymean - (b1 * xmean)

        # Estimated y values
        y_est = [(b1 * x) + b_nought for x in x_var]

        # r-squared numerator
        y_est_sum = sum((y - ymean) ** 2 for y in y_est)

        # Standard error
        y_err_sum = sum((y_est[i] - y_var[i]) ** 2 for i in range(len(y_var)))
        std_err = (y_err_sum / (len(y_var) - 2)) ** 0.5

        # Residuals
        residuals = [round(y_var[i] - y_est[i], decimals) for i in range(len(y_var))]

        # Return dictionary
        return {
            "intercept": round(b_nought, decimals),
            "slope": round(b1, decimals),
            "r_squared": round(y_est_sum / y_sqr_sum, decimals),
            "std_error": round(std_err, decimals),
            "residuals": residuals
        }

    def linear_regression_predict(self, x_values, model, decimals=2):
        """Generate predicted y-values using a linear regression model.

        Parameters:

        x_values (list): the list of independent variable values to generate predictions for. The input list will not be altered by the operation.

        model (dict): a dictionary returned by the linear_regression function containing at least the keys "intercept" and "slope".

        decimals (integer): the number of decimal places predicted values will be rounded to."""
        # Validate x_values
        if not isinstance(x_values, list) or not x_values:
            raise TypeError("'x_values' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in x_values):
            raise ValueError("All elements in 'x_values' must be numeric.")

        # Validate model structure
        if not isinstance(model, dict):
            raise TypeError("'model' must be a dictionary returned by linear_regression.")
        if "intercept" not in model or "slope" not in model:
            raise ValueError("Model dictionary must contain 'intercept' and 'slope' keys.")

        b0 = model["intercept"]
        b1 = model["slope"]

        return [round(b0 + b1 * x, decimals) for x in x_values]

    def multiple_linear_regression(self, x_vars, y_var, decimals=2):
        """Perform multiple linear regression using a list of predictor variables ('x_vars') and a dependent variable ('y_var'). Returns a dictionary containing the intercept (key: intercept), coefficients (key: coefficients), slope (key: slope), r-squared (key: r_squared), standard error (key: std_error), and residuals (key: residuals), as well as adjusted r-squared (key: adjusted_r_squared), an ANOVA table (key: anova), and coefficient-level t-test results (key: coefficients_significance).

        Parameters:

        x_vars (list): a nested list where each sub-list represents one independent variable the operation will be performed on. All predictor lists must be the same length. The input list will not be altered by the operation

        y_var (list): the list of values representing the dependent variable the operation will be performed on. The input list will not be altered by the operation.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(x_vars, list) or not x_vars:
            raise TypeError("'X_vars' must be a non-empty list of predictor lists.")
        if not all(isinstance(col, list) and col for col in x_vars):
            raise TypeError("Each predictor must be a non-empty list.")
        if not isinstance(y_var, list) or not y_var:
            raise TypeError("'y_var' must be a non-empty list.")

        n = len(y_var)
        k = len(x_vars)

        # All predictors must match y_var length
        for col in x_vars:
            if len(col) != n:
                raise ValueError("All predictor lists must match the length of 'y_var'.")
            if not all(isinstance(v, (int, float)) for v in col):
                raise ValueError("All predictor values must be numeric.")

        if not all(isinstance(v, (int, float)) for v in y_var):
            raise ValueError("All elements in 'y_var' must be numeric.")

        # Build design matrix X with intercept column
        X = []
        for i in range(n):
            row = [1]  # intercept
            for col in x_vars:
                row.append(col[i])
            X.append(row)

        # Compute X^T X
        XT = list(zip(*X))
        XTX = [[sum(a*b for a, b in zip(row_i, col_j))
                for col_j in zip(*X)] for row_i in XT]

        #  Compute X^T y 
        XTy = [sum(x_i * y_var[i] for i, x_i in enumerate(col)) for col in XT]

        # Invert X^T X (Gauss-Jordan elimination)
        size = len(XTX)
        aug = [XTX[i] + [XTy[i]] for i in range(size)]

        for i in range(size):
            pivot = aug[i][i]
            if pivot == 0:
                raise ValueError("Matrix inversion failed; predictors may be collinear.")
            for j in range(i, size + 1):
                aug[i][j] /= pivot
            for r in range(size):
                if r != i:
                    factor = aug[r][i]
                    for c in range(i, size + 1):
                        aug[r][c] -= factor * aug[i][c]

        # Extract coefficients
        beta = [aug[i][-1] for i in range(size)]

        b0 = round(beta[0], decimals)
        coeffs = [round(b, decimals) for b in beta[1:]]

        # Predictions
        predicted = []
        for i in range(n):
            y_hat = beta[0]
            for j in range(k):
                y_hat += beta[j+1] * x_vars[j][i]
            predicted.append(round(y_hat, decimals))
        
        # Residuals
        residuals = [round(y_var[i] - predicted[i], decimals) for i in range(n)]

        # R-squared
        y_mean = sum(y_var) / n
        ss_tot = sum((y - y_mean)**2 for y in y_var)
        ss_res = sum((y_var[i] - predicted[i])**2 for i in range(n))
        r_squared = round(1 - (ss_res / ss_tot), decimals)

        # Adjusted R-squared
        if n - k - 1 > 0:
            adj_r_squared = 1 - ((1 - r_squared) * (n - 1) / (n - k - 1))
            adj_r_squared = round(adj_r_squared, decimals)
        else:
            adj_r_squared = None  # undefined when n <= k + 1

        # Standard error of regression
        std_error = round((ss_res / (n - k - 1)) ** 0.5, decimals)

        # Regression ANOVA
        ss_reg = ss_tot - ss_res
        df_reg = k
        df_res = n - k - 1
        df_tot = n - 1

        ms_reg = ss_reg / df_reg
        ms_res = ss_res / df_res

        F_stat = ms_reg / ms_res
        F_stat = round(F_stat, decimals)

        # Significance using F-table
        anova_significance = self._check_significance(
            stat_value=F_stat,
            df=(df_reg, df_res),
            test_type='f',
            alpha=0.05,
            tail_type='two_tail'
        )

        anova_table = {
            "ss_reg": round(ss_reg, decimals),
            "ss_res": round(ss_res, decimals),
            "ss_tot": round(ss_tot, decimals),
            "df_reg": df_reg,
            "df_res": df_res,
            "df_tot": df_tot,
            "ms_reg": round(ms_reg, decimals),
            "ms_res": round(ms_res, decimals),
            "F_stat": F_stat,
            "significance": anova_significance
        }

        # Coefficient-level t-tests
        # Standard errors for coefficients come from diagonal of (X^T X)^(-1) * MSE
        mse = ms_res  # same as residual MS

        # Extract inverse(X^T X) from augmented matrix
        inv_xtx = [[aug[i][j] for j in range(size)] for i in range(size)]

        coeffs_significance = {}
        for idx in range(1, size):  # skip intercept at index 0
            var_name = f"X{idx}"
            coef_value = beta[idx]
            coef_se = (mse * inv_xtx[idx][idx]) ** 0.5
            t_stat = coef_value / coef_se if coef_se != 0 else float('inf')

            t_stat_rounded = round(t_stat, decimals)
            coef_se_rounded = round(coef_se, decimals)

            significance = self._check_significance(
                stat_value=t_stat_rounded,
                df=df_res,
                test_type='t',
                alpha=0.05,
                tail_type='two_tail'
            )

            coeffs_significance[var_name] = {
                "coefficient": round(coef_value, decimals),
                "std_error": coef_se_rounded,
                "t_stat": t_stat_rounded,
                "significance": significance
            }

        return {
            "intercept": b0,
            "coefficients": coeffs,
            "predicted": predicted,
            "residuals": residuals,
            "r_squared": r_squared,
            "adjusted_r_squared": adj_r_squared,
            "std_error": std_error,
            "anova": anova_table,
            "coefficients_significance": coeffs_significance
        }
    
    def multiple_regression_predict(self, x_vars, model, decimals=2):
        """Generate predicted y-values using a multiple linear regression model.

        Parameters:

        x_vars (list): a nested list where each sub-list represents one independent variable to generate predictions for. All predictor lists must be the same length. The input list will not be altered by the operation.

        model (dict): a dictionary returned by the multiple_linear_regression function containing the keys 'intercept' and 'coefficients'.

        decimals (integer): number of decimal places to round returned values to."""
        # Validate x_vars
        if not isinstance(x_vars, list) or not x_vars:
            raise TypeError("'x_vars' must be a non-empty list of predictor lists.")
        if not all(isinstance(col, list) and col for col in x_vars):
            raise TypeError("Each predictor in 'x_vars' must be a non-empty list.")

        n = len(x_vars[0])
        for col in x_vars:
            if len(col) != n:
                raise ValueError("All predictor lists in 'x_vars' must be the same length.")
            if not all(isinstance(v, (int, float)) for v in col):
                raise ValueError("All predictor values must be numeric.")

        # Validate model
        if not isinstance(model, dict):
            raise TypeError("'model' must be a dictionary returned by multiple_linear_regression.")
        if "intercept" not in model or "coefficients" not in model:
            raise ValueError("Model dictionary must contain 'intercept' and 'coefficients' keys.")

        b0 = model["intercept"]
        coeffs = model["coefficients"]

        if len(coeffs) != len(x_vars):
            raise ValueError("Number of predictors in 'x_vars' must match the model's coefficient count.")

        # Generate predictions
        predicted = []
        for i in range(n):
            y_hat = b0
            for j in range(len(coeffs)):
                y_hat += coeffs[j] * x_vars[j][i]
            predicted.append(round(y_hat, decimals))

        return predicted

    def ci_regression_slope(self, slope, std_error, df, alpha=0.05, decimals=2):
        """Calculate a confidence interval for a regression slope ('slope') using its standard error ('std_error') and degrees of freedom ('df'). Returns a dictionary with the following keys: slope, lower, upper,  and margin_of_error.

        Parameters:

        slope (float): the estimated regression slope coefficient.

        std_error (float): the standard error associated with the slope estimate.

        df (integer or float): the degrees of freedom for the regression model.

        alpha (float): the significance level for the confidence interval.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(slope, (int, float)):
            raise TypeError("'slope' must be a numeric value.")
        if not isinstance(std_error, (int, float)) or std_error < 0:
            raise TypeError("'std_error' must be a non-negative numeric value.")
        if not isinstance(df, (int, float)) or df <= 0:
            raise TypeError("'df' must be a positive numeric value.")
        if not isinstance(alpha, float):
            raise TypeError("'alpha' must be a float value.")

        # t-critical value (two-tailed)
        t_crit = self._get_critical_value(
            stat_type='t',
            alpha=alpha,
            df=df
        )

        # Margin of error
        moe = t_crit * std_error

        lower = slope - moe
        upper = slope + moe

        return {
            "slope": round(slope, decimals),
            "lower": round(lower, decimals),
            "upper": round(upper, decimals),
            "margin_of_error": round(moe, decimals)
        }

    def vif(self, predictors, predictor_names=None, decimals=2):
        """Calculate Variance Inflation Factors (VIF) for a set of predictor variables. Returns a dictionary mapping each predictor to its VIF value.

        Parameters:

        predictors (list): a nested list where each element is a list of numeric values representing one predictor variable. All predictors must be the same length.

        predictor_names (list): optional list of names corresponding to each predictor. If None, predictors will be labeled by index.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(predictors, list) or len(predictors) < 2:
            raise TypeError("'predictors' must be a list containing at least two predictors.")
        length_set = {len(p) for p in predictors}
        if len(length_set) != 1:
            raise ValueError("All predictors must have the same number of observations.")
        for p in predictors:
            if not all(isinstance(v, (int, float)) for v in p):
                raise ValueError("All predictor values must be numeric.")

        k = len(predictors)

        # Handle predictor names
        if predictor_names is None:
            predictor_names = [f"X{i+1}" for i in range(k)]
        if len(predictor_names) != k:
            raise ValueError("'predictor_names' must match the number of predictors.")

        vif_values = {}

        # Compute VIF for each predictor
        for i in range(k):
            # Dependent variable = predictor i
            y = predictors[i]

            # Independent variables = all other predictors
            X = [p for j, p in enumerate(predictors) if j != i]

            # Run regression using your existing engine
            result = self.multiple_linear_regression(X, y)

            # Extract R^2
            r2 = result["r_squared"]

            # Compute VIF
            if r2 >= 1:
                vif = float("inf")
            else:
                vif = 1 / (1 - r2)

            vif_values[predictor_names[i]] = round(vif, decimals)

        return vif_values

    def lowess(self, x_vals, y_vals, fraction=0.3, decimals=2):
        """Perform LOWESS (Locally Weighted Scatterplot Smoothing) on paired data ('x_vals', 'y_vals'). Returns a dictionary with the following keys: x and y_smooth.

        Parameters:

        x_vals (list): the list of x-values. The input list will not be altered by the operation.

        y_vals (list): the list of y-values corresponding to x_vals. The input list will not be altered.

        fraction (float): the fraction of the data used in each local regression (smoothing parameter).

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(x_vals, list) or not isinstance(y_vals, list):
            raise TypeError("'x_vals' and 'y_vals' must be lists.")
        if len(x_vals) != len(y_vals) or len(x_vals) < 3:
            raise ValueError("'x_vals' and 'y_vals' must be the same length and contain at least 3 points.")
        if not all(isinstance(v, (int, float)) for v in x_vals):
            raise ValueError("'x_vals' must contain only numeric values.")
        if not all(isinstance(v, (int, float)) for v in y_vals):
            raise ValueError("'y_vals' must contain only numeric values.")
        if not isinstance(fraction, float) or not (0 < fraction <= 1):
            raise ValueError("'fraction' must be a float between 0 and 1.")

        n = len(x_vals)
        k = max(2, int(fraction * n))  # neighborhood size

        # Sort by x
        sorted_pairs = sorted(zip(x_vals, y_vals), key=lambda p: p[0])
        xs_sorted = [p[0] for p in sorted_pairs]
        ys_sorted = [p[1] for p in sorted_pairs]

        y_smooth_sorted = []

        for i in range(n):
            x0 = xs_sorted[i]

            # Compute distances
            distances = [abs(x - x0) for x in xs_sorted]

            # Identify k nearest neighbors
            idx_sorted = sorted(range(n), key=lambda j: distances[j])
            neighbors = idx_sorted[:k]

            # Maximum distance in neighborhood
            d_max = distances[neighbors[-1]] if distances[neighbors[-1]] != 0 else 1e-12

            # Compute tricube weights
            weights = []
            for j in neighbors:
                u = distances[j] / d_max
                w = (1 - u**3)**3 if u < 1 else 0
                weights.append(w)

            # Weighted linear regression (simple, explicit)
            X = [xs_sorted[j] for j in neighbors]
            Y = [ys_sorted[j] for j in neighbors]
            W = weights

            # Weighted means
            w_sum = sum(W)
            x_bar = sum(w * x for w, x in zip(W, X)) / w_sum
            y_bar = sum(w * y for w, y in zip(W, Y)) / w_sum

            # Weighted slope
            num = sum(w * (x - x_bar) * (y - y_bar) for w, x, y in zip(W, X, Y))
            den = sum(w * (x - x_bar)**2 for w, x in zip(W, X))
            slope = num / den if den != 0 else 0

            # Weighted intercept
            intercept = y_bar - slope * x_bar

            # Smoothed prediction
            y0 = intercept + slope * x0
            y_smooth_sorted.append(y0)

        # Reorder to match original x order
        mapping = {xs_sorted[i]: y_smooth_sorted[i] for i in range(n)}
        y_smooth_original = [mapping[x] for x in x_vals]

        return {
            "x": x_vals,
            "y_smooth": [round(v, decimals) for v in y_smooth_original]
        }

    def loess(self, x_vals, y_vals, fraction=0.3, degree=2, alpha=0.05, decimals=2):
        """Perform LOESS (Locally Estimated Scatterplot Smoothing) on paired data ('x_vals', 'y_vals'), returning smoothed values and confidence bands. Returns a dictionary with the following keys: x, y_smooth, standard_error, lower, and upper.

        Parameters:

        x_vals (list): the list of x-values. The input list will not be altered by the operation.

        y_vals (list): the list of y-values corresponding to x_vals. The input list will not be altered by the operation.

        fraction (float): the fraction of the data used in each local regression (smoothing parameter).

        degree (integer): the degree of the local polynomial (0, 1, or 2). Default is 2 (quadratic LOESS).

        alpha (float): the significance level for the confidence interval.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(x_vals, list) or not isinstance(y_vals, list):
            raise TypeError("'x_vals' and 'y_vals' must be lists.")
        if len(x_vals) != len(y_vals) or len(x_vals) < 3:
            raise ValueError("'x_vals' and 'y_vals' must be the same length and contain at least 3 points.")
        if not all(isinstance(v, (int, float)) for v in x_vals):
            raise ValueError("'x_vals' must contain only numeric values.")
        if not all(isinstance(v, (int, float)) for v in y_vals):
            raise ValueError("'y_vals' must contain only numeric values.")
        if not isinstance(fraction, float) or not (0 < fraction <= 1):
            raise ValueError("'fraction' must be a float between 0 and 1.")
        if degree not in (0, 1, 2):
            raise ValueError("'degree' must be 0, 1, or 2.")
        if not isinstance(alpha, float):
            raise TypeError("'alpha' must be a float value.")

        n = len(x_vals)
        k = max(2, int(fraction * n))  # neighborhood size

        # Sort by x
        sorted_pairs = sorted(zip(x_vals, y_vals), key=lambda p: p[0])
        xs_sorted = [p[0] for p in sorted_pairs]
        ys_sorted = [p[1] for p in sorted_pairs]

        y_smooth_sorted = []
        se_sorted = []

        p = degree + 1  # number of coefficients

        for i in range(n):
            x0 = xs_sorted[i]

            # Compute distances
            distances = [abs(x - x0) for x in xs_sorted]

            # Identify k nearest neighbors
            idx_sorted = sorted(range(n), key=lambda j: distances[j])
            neighbors = idx_sorted[:k]

            # Maximum distance in neighborhood
            d_max = distances[neighbors[-1]] if distances[neighbors[-1]] != 0 else 1e-12

            # Compute tricube weights
            weights = []
            for j in neighbors:
                u = distances[j] / d_max
                w = (1 - u**3)**3 if u < 1 else 0
                weights.append(w)

            # Build weighted design matrix for polynomial regression
            X = []
            Y = []
            W = weights

            for j in neighbors:
                x = xs_sorted[j]
                y = ys_sorted[j]
                if degree == 0:
                    X.append([1])
                elif degree == 1:
                    X.append([1, x])
                else:  # degree == 2
                    X.append([1, x, x**2])
                Y.append(y)

            # Build normal equations: (X^T W X) and (X^T W y)
            XT_W_X = [[0]*p for _ in range(p)]
            XT_W_y = [0]*p

            for row, w, y in zip(X, W, Y):
                for a in range(p):
                    XT_W_y[a] += w * row[a] * y
                    for b in range(p):
                        XT_W_X[a][b] += w * row[a] * row[b]

            # Solve (X^T W X) b = X^T W y using Gaussian elimination
            A = [XT_W_X[i] + [XT_W_y[i]] for i in range(p)]

            # Forward elimination
            for col in range(p):
                pivot = A[col][col]
                if abs(pivot) < 1e-12:
                    continue
                for row in range(col+1, p):
                    factor = A[row][col] / pivot
                    for c in range(col, p+1):
                        A[row][c] -= factor * A[col][c]

            # Back substitution
            coeffs = [0]*p
            for row in range(p-1, -1, -1):
                rhs = A[row][p]
                for c in range(row+1, p):
                    rhs -= A[row][c] * coeffs[c]
                coeffs[row] = rhs / (A[row][row] if abs(A[row][row]) > 1e-12 else 1e-12)

            # Predict smoothed value
            if degree == 0:
                y0 = coeffs[0]
                h = [1]
            elif degree == 1:
                y0 = coeffs[0] + coeffs[1] * x0
                h = [1, x0]
            else:
                y0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0**2
                h = [1, x0, x0**2]

            y_smooth_sorted.append(y0)

            # Compute local residual variance
            residuals = []
            for row, y, w in zip(X, Y, W):
                y_hat = sum(c * r for c, r in zip(coeffs, row))
                residuals.append(w * (y - y_hat)**2)

            w_sum = sum(W)
            df_local = max(1, w_sum - p)
            sigma2 = sum(residuals) / df_local

            # Invert (X^T W X) for variance calculation
            # Build augmented matrix for inversion
            M = [XT_W_X[i] + [1 if i == j else 0 for j in range(p)] for i in range(p)]

            # Gaussian elimination for matrix inversion
            # Forward elimination
            for col in range(p):
                pivot = M[col][col]
                if abs(pivot) < 1e-12:
                    continue
                for row in range(col+1, p):
                    factor = M[row][col] / pivot
                    for c in range(col, 2*p):
                        M[row][c] -= factor * M[col][c]

            # Back substitution
            for col in range(p-1, -1, -1):
                pivot = M[col][col]
                if abs(pivot) < 1e-12:
                    continue
                for row in range(col-1, -1, -1):
                    factor = M[row][col] / pivot
                    for c in range(col, 2*p):
                        M[row][c] -= factor * M[col][c]

            # Normalize rows
            for i in range(p):
                pivot = M[i][i]
                if abs(pivot) < 1e-12:
                    continue
                for c in range(i, 2*p):
                    M[i][c] /= pivot

            # Extract inverse matrix
            XT_W_X_inv = [row[p:] for row in M]

            # Variance of prediction: h^T (X^T W X)^(-1) h * sigma^2
            var_pred = 0
            for a in range(p):
                for b in range(p):
                    var_pred += h[a] * XT_W_X_inv[a][b] * h[b]
            var_pred *= sigma2

            se = var_pred**0.5
            se_sorted.append(se)

        # t-critical value
        t_crit = self._get_critical_value(stat_type='t', alpha=alpha, df=n - p)

        # Confidence bands
        lower_sorted = [y - t_crit * s for y, s in zip(y_smooth_sorted, se_sorted)]
        upper_sorted = [y + t_crit * s for y, s in zip(y_smooth_sorted, se_sorted)]

        # Reorder to match original x order
        mapping = {
            xs_sorted[i]: (
                y_smooth_sorted[i],
                se_sorted[i],
                lower_sorted[i],
                upper_sorted[i]
            )
            for i in range(n)
        }

        y_smooth_original = []
        se_original = []
        lower_original = []
        upper_original = []

        for x in x_vals:
            y_s, se, lo, up = mapping[x]
            y_smooth_original.append(round(y_s, decimals))
            se_original.append(round(se, decimals))
            lower_original.append(round(lo, decimals))
            upper_original.append(round(up, decimals))

        return {
            "x": x_vals,
            "y_smooth": y_smooth_original,
            "standard_error": se_original,
            "lower": lower_original,
            "upper": upper_original
        }

    # ============================================================
    # HYPOTHESIS TESTING - MEANS
    # ============================================================

    def one_sample_t_test(self, vals, mu, alpha=0.05, decimals=2):
        """Calculate the one-sample T value for an input list ('val's) to test for a difference of the mean from a theoretical mean ('mu'). Returns a dictionary containing the t-statistic (key: t_stat), degrees of freedom (key: df), significance decision (key: significance), and the critical value used (key: critical_value).

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        mu (integer,float): the theoretical mean.

        alpha (float): the significance level for the t-test.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        if not isinstance(mu, (int, float)):
            raise TypeError("'mu' must be an integer or float.")
        if not isinstance(alpha, float):
            raise TypeError("'alpha' must be a float value.")

        n = len(vals)
        df = n - 1

        # Mean
        mean = sum(vals) / n

        # Standard deviation (sample)
        sum_sqr = 0
        for v in vals:
            sum_sqr += (v - mean) ** 2
        std_dev = math.sqrt(sum_sqr / df)

        # Correct t-statistic formula
        t_stat = (mean - mu) / (std_dev / math.sqrt(n))
        t_stat = round(t_stat, decimals)

        # Significance check (two-tailed by default)
        significance_result = self._check_significance(
            stat_value=t_stat,
            df=df,
            test_type='t',
            alpha=alpha,
            tail_type='two_tail'
        )

        # Extract critical value and df_used from the returned string
        try:
            crit_part = significance_result.split("Critical Value:")[1]
            critical_value = float(crit_part.split("(")[0].strip())
            df_used = crit_part.split("DF used:")[1].replace(").", "").strip()
        except Exception:
            critical_value = None
            df_used = None

        return {
            "t_stat": t_stat,
            "df": df,
            "alpha": alpha,
            "critical_value": critical_value,
            "df_used": df_used,
            "significance": significance_result
        }

    def unpaired_t_test(self, x_var, y_var, alpha=0.05, decimals=2):
        """Calculate the unpaired (independent samples) T value for two variables ('x_var', 'y_var') to test for a difference of means. Returns a dictionary containing the t-statistic (key: t_stat), degrees of freedom (key: df), significance decision (key: significance), and the critical value used (key: critical_value).

        Parameters:

        x_var (list): the list of values representing the x variable the operation will be performed on. The input list will not be altered by the operation.

        y_var (list): the list of values representing the y variable the operation will be performed on. The input list will not be altered by the operation.

        alpha (float): the significance level for the t-test.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(x_var, list) or not x_var:
            raise TypeError("'x_var' must be a non-empty list.")
        if not isinstance(y_var, list) or not y_var:
            raise TypeError("'y_var' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in x_var):
            raise ValueError("All elements in 'x_var' must be numeric.")
        if not all(isinstance(v, (int, float)) for v in y_var):
            raise ValueError("All elements in 'y_var' must be numeric.")
        if not isinstance(alpha, float):
            raise TypeError("'alpha' must be a float value.")

        nx = len(x_var)
        ny = len(y_var)
        df = nx + ny - 2

        # Means
        xmean = sum(x_var) / nx
        ymean = sum(y_var) / ny

        # Sum of squared deviations
        sum_sqr_x = sum((v - xmean) ** 2 for v in x_var)
        sum_sqr_y = sum((v - ymean) ** 2 for v in y_var)

        # Correct pooled variance
        s2 = (sum_sqr_x + sum_sqr_y) / df

        # T-statistic
        t_stat = (xmean - ymean) / math.sqrt(s2 * (1/nx + 1/ny))
        t_stat = round(t_stat, decimals)

        # Significance check (two-tailed)
        significance_result = self._check_significance(
            stat_value=t_stat,
            df=df,
            test_type='t',
            alpha=alpha,
            tail_type='two_tail'
        )

        # Extract critical value and df_used
        try:
            crit_part = significance_result.split("Critical Value:")[1]
            critical_value = float(crit_part.split("(")[0].strip())
            df_used = crit_part.split("DF used:")[1].replace(").", "").strip()
        except Exception:
            critical_value = None
            df_used = None

        return {
            "t_stat": t_stat,
            "df": df,
            "alpha": alpha,
            "critical_value": critical_value,
            "df_used": df_used,
            "significance": significance_result
        }

    def paired_t_test(self, x_var, y_var, alpha=0.05, decimals=2):
        """Calculate the paired-samples T value for two dependent variables ('x_var', 'y_var') to test for a difference of means. Returns a dictionary containing the t-statistic (key: t_stat), degrees of freedom (key: df), significance decision (key: significance), and the critical value used (key: critical_value).

        Parameters:

        x_var (list): the list of values representing the x variable the operation will be performed on. The input list will not be altered by the operation.

        y_var (list): the list of values representing the y variable the operation will be performed on. The input list will not be altered by the operation.

        alpha (float): the significance level for the t-test.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(x_var, list) or not x_var:
            raise TypeError("'x_var' must be a non-empty list.")
        if not isinstance(y_var, list) or not y_var:
            raise TypeError("'y_var' must be a non-empty list.")
        if len(x_var) != len(y_var):
            raise ValueError("'x_var' and 'y_var' must be the same length.")
        if not all(isinstance(v, (int, float)) for v in x_var):
            raise ValueError("All elements in 'x_var' must be numeric.")
        if not all(isinstance(v, (int, float)) for v in y_var):
            raise ValueError("All elements in 'y_var' must be numeric.")
        if not isinstance(alpha, float):
            raise TypeError("'alpha' must be a float value.")

        # Differences
        d = [x_var[i] - y_var[i] for i in range(len(x_var))]
        n = len(d)
        df = n - 1

        # Mean of differences
        mean_d = sum(d) / n

        # Standard deviation of differences (sample)
        sum_sqr = sum((val - mean_d) ** 2 for val in d)
        std_dev = math.sqrt(sum_sqr / df)

        # T-statistic
        t_stat = mean_d / (std_dev / math.sqrt(n))
        t_stat = round(t_stat, decimals)

        # Significance check (two-tailed)
        significance_result = self._check_significance(
            stat_value=t_stat,
            df=df,
            test_type='t',
            alpha=alpha,
            tail_type='two_tail'
        )

        # Extract critical value and df_used
        try:
            crit_part = significance_result.split("Critical Value:")[1]
            critical_value = float(crit_part.split("(")[0].strip())
            df_used = crit_part.split("DF used:")[1].replace(").", "").strip()
        except Exception:
            critical_value = None
            df_used = None

        return {
            "t_stat": t_stat,
            "df": df,
            "alpha": alpha,
            "critical_value": critical_value,
            "df_used": df_used,
            "significance": significance_result
        }

    def welch_t_test(self, x_var, y_var, alpha=0.05, decimals=2):
        """Calculate Welch's two-sample T value for two independent variables (x_var, y_var) without assuming equal variances. Returns a dictionary containing the t-statistic (key: t_stat), degrees of freedom (key: df), significance decision (key: significance), and the critical value used (key: critical_value).

        Parameters:

        x_var (list): the list of values representing the x variable the operation will be performed on. The input list will not be altered by the operation.

        y_var (list): the list of values representing the y variable the operation will be performed on. The input list will not be altered by the operation.

        alpha (float): the significance level for the t-test.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(x_var, list) or not x_var:
            raise TypeError("'x_var' must be a non-empty list.")
        if not isinstance(y_var, list) or not y_var:
            raise TypeError("'y_var' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in x_var):
            raise ValueError("All elements in 'x_var' must be numeric.")
        if not all(isinstance(v, (int, float)) for v in y_var):
            raise ValueError("All elements in 'y_var' must be numeric.")
        if not isinstance(alpha, float):
            raise TypeError("'alpha' must be a float value.")

        nx = len(x_var)
        ny = len(y_var)

        # Means
        xmean = sum(x_var) / nx
        ymean = sum(y_var) / ny

        # Sample variances
        var_x = sum((v - xmean) ** 2 for v in x_var) / (nx - 1)
        var_y = sum((v - ymean) ** 2 for v in y_var) / (ny - 1)

        # Welch t-statistic
        t_stat = (xmean - ymean) / math.sqrt(var_x / nx + var_y / ny)
        t_stat = round(t_stat, decimals)

        # Welch–Satterthwaite degrees of freedom
        numerator = (var_x / nx + var_y / ny) ** 2
        denominator = ((var_x / nx) ** 2) / (nx - 1) + ((var_y / ny) ** 2) / (ny - 1)
        df = numerator / denominator
        df = int(df)  # Conservative rounding down for table lookup

        # Significance check (two-tailed)
        significance_result = self._check_significance(
            stat_value=t_stat,
            df=df,
            test_type='t',
            alpha=alpha,
            tail_type='two_tail'
        )

        # Extract critical value and df_used
        try:
            crit_part = significance_result.split("Critical Value:")[1]
            critical_value = float(crit_part.split("(")[0].strip())
            df_used = crit_part.split("DF used:")[1].replace(").", "").strip()
        except Exception:
            critical_value = None
            df_used = None

        return {
            "t_stat": t_stat,
            "df": df,
            "alpha": alpha,
            "critical_value": critical_value,
            "df_used": df_used,
            "significance": significance_result
        }

    def ci_mean(self, vals, alpha=0.05, decimals=2):
        """Calculate a confidence interval for the mean of a numeric sample ('vals'). Returns a dictionary with the following keys: mean, lower, upper, and margin_of_error.

        Parameters:

        vals (list): the list of numeric values representing the sample. The input list will not be altered by the operation.

        alpha (float): the significance level for the confidence interval.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("'vals' must contain only numeric values.")
        if not isinstance(alpha, float):
            raise TypeError("'alpha' must be a float value.")

        n = len(vals)
        mean_val = sum(vals) / n

        # Sample standard deviation
        var = sum((x - mean_val)**2 for x in vals) / (n - 1)
        sd = var ** 0.5

        # Standard error
        se = sd / (n ** 0.5)

        # Degrees of freedom
        df = n - 1

        # t-critical value (two-tailed)
        t_crit = self._get_critical_value(
            stat_type='t',
            alpha=alpha,
            df=df
        )

        # Margin of error
        moe = t_crit * se

        lower = mean_val - moe
        upper = mean_val + moe

        return {
            "mean": round(mean_val, decimals),
            "lower": round(lower, decimals),
            "upper": round(upper, decimals),
            "margin_of_error": round(moe, decimals)
        }
 
    def ci_difference_of_means(self, vals1, vals2, alpha=0.05, decimals=2):
        """Calculate a confidence interval for the difference of means between two independent samples ('vals1' and 'vals2') using Welch's method. Returns a dictionary with the following keys:  difference, lower, upper, and margin_of_error.

        Parameters:

        vals1 (list): the list of numeric values representing the first sample. The input list will not be altered by the operation.

        vals2 (list): the list of numeric values representing the second sample. The input list will not be altered by the operation.

        alpha (float): the significance level for the confidence interval.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(vals1, list) or not vals1:
            raise TypeError("'vals1' must be a non-empty list.")
        if not isinstance(vals2, list) or not vals2:
            raise TypeError("'vals2' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals1):
            raise ValueError("'vals1' must contain only numeric values.")
        if not all(isinstance(v, (int, float)) for v in vals2):
            raise ValueError("'vals2' must contain only numeric values.")
        if not isinstance(alpha, float):
            raise TypeError("'alpha' must be a float value.")

        n1 = len(vals1)
        n2 = len(vals2)

        mean1 = sum(vals1) / n1
        mean2 = sum(vals2) / n2

        # Sample variances
        var1 = sum((x - mean1)**2 for x in vals1) / (n1 - 1)
        var2 = sum((x - mean2)**2 for x in vals2) / (n2 - 1)

        # Standard error of the difference
        se = ((var1 / n1) + (var2 / n2)) ** 0.5

        # Welch–Satterthwaite degrees of freedom
        df_num = (var1 / n1 + var2 / n2) ** 2
        df_den = ((var1 / n1) ** 2) / (n1 - 1) + ((var2 / n2) ** 2) / (n2 - 1)
        df = df_num / df_den if df_den else 0

        # t-critical value (two-tailed)
        t_crit = self._get_critical_value(
            stat_type='t',
            alpha=alpha,
            df=df
        )

        # Difference of means
        diff = mean1 - mean2

        # Margin of error
        moe = t_crit * se

        lower = diff - moe
        upper = diff + moe

        return {
            "difference": round(diff, decimals),
            "lower": round(lower, decimals),
            "upper": round(upper, decimals),
            "margin_of_error": round(moe, decimals)
        }

    # ============================================================
    # HYPOTHESIS TESTING - PROPORTIONS
    # ============================================================

    def one_proportion_z_test(self, successes, n, p0, alpha=0.05, decimals=2):
        """Perform a one-proportion Z test to compare an observed proportion (successes / n) against a theoretical proportion (p0). Returns a dictionary containing the z-statistic (key: z_stat), degrees of freedom used for lookup (key: df_used), significance decision (key: significance), and the critical value used (key: critical_value).

        Parameters:

        successes (integer): the number of observed successes.

        n (integer): the total number of trials.

        p0 (float): the theoretical population proportion.

        alpha (float): the significance level for the z-test.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(successes, int) or successes < 0:
            raise TypeError("'successes' must be a non-negative integer.")
        if not isinstance(n, int) or n <= 0:
            raise TypeError("'n' must be a positive integer.")
        if successes > n:
            raise ValueError("'successes' cannot exceed 'n'.")
        if not isinstance(p0, float):
            raise TypeError("'p0' must be a float value.")
        if not isinstance(alpha, float):
            raise TypeError("'alpha' must be a float value.")

        # Observed proportion
        p_hat = successes / n

        # Standard error
        se = math.sqrt(p0 * (1 - p0) / n)

        # Z-statistic
        z_stat = (p_hat - p0) / se
        z_stat = round(z_stat, decimals)

        # Use Z-distribution row ('max') in T-table
        significance_result = self._check_significance(
            stat_value=z_stat,
            df=999999,  # triggers 'max' row
            test_type='t',
            alpha=alpha,
            tail_type='two_tail'
        )

        # Extract critical value and df_used
        try:
            crit_part = significance_result.split("Critical Value:")[1]
            critical_value = float(crit_part.split("(")[0].strip())
            df_used = crit_part.split("DF used:")[1].replace(").", "").strip()
        except Exception:
            critical_value = None
            df_used = None

        return {
            "z_stat": z_stat,
            "n": n,
            "successes": successes,
            "p_hat": round(p_hat, decimals),
            "p0": p0,
            "alpha": alpha,
            "critical_value": critical_value,
            "df_used": df_used,
            "significance": significance_result
        }

    def two_proportion_z_test(self, successes1, n1, successes2, n2, alpha=0.05, decimals=2):
        """Perform a two-proportion Z test to compare the proportions of two independent samples. Returns a dictionary containing the z-statistic (key: z_stat), degrees of freedom used for lookup (key: df_used), significance decision (key: significance), and the critical value used (key: critical_value).

        Parameters:

        successes1 (integer): number of successes in sample 1.

        n1 (integer): total trials in sample 1.

        successes2 (integer): number of successes in sample 2.

        n2 (integer): total trials in sample 2.

        alpha (float): the significance level for the z-test.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        for name, val in [('successes1', successes1), ('successes2', successes2)]:
            if not isinstance(val, int) or val < 0:
                raise TypeError(f"'{name}' must be a non-negative integer.")

        for name, val in [('n1', n1), ('n2', n2)]:
            if not isinstance(val, int) or val <= 0:
                raise TypeError(f"'{name}' must be a positive integer.")

        if successes1 > n1 or successes2 > n2:
            raise ValueError("Success counts cannot exceed sample sizes.")

        if not isinstance(alpha, float):
            raise TypeError("'alpha' must be a float value.")

        # Sample proportions
        p1 = successes1 / n1
        p2 = successes2 / n2

        # Pooled proportion
        pooled = (successes1 + successes2) / (n1 + n2)

        # Standard error
        se = math.sqrt(pooled * (1 - pooled) * (1/n1 + 1/n2))

        # Z-statistic
        z_stat = (p1 - p2) / se
        z_stat = round(z_stat, decimals)

        # Use Z-distribution row ('max') in T-table
        significance_result = self._check_significance(
            stat_value=z_stat,
            df=999999,  # triggers 'max' row
            test_type='t',
            alpha=alpha,
            tail_type='two_tail'
        )

        # Extract critical value and df_used
        try:
            crit_part = significance_result.split("Critical Value:")[1]
            critical_value = float(crit_part.split("(")[0].strip())
            df_used = crit_part.split("DF used:")[1].replace(").", "").strip()
        except Exception:
            critical_value = None
            df_used = None

        return {
            "z_stat": z_stat,
            "p1": round(p1, decimals),
            "p2": round(p2, decimals),
            "pooled": round(pooled, decimals),
            "n1": n1,
            "n2": n2,
            "alpha": alpha,
            "critical_value": critical_value,
            "df_used": df_used,
            "significance": significance_result
        }

    def ci_proportion(self, vals, alpha=0.05, decimals=2):
        """Calculate a confidence interval for a proportion based on a binary sample ('vals'). Returns a dictionary with the following keys: proportion, lower, upper, and margin_of_error.

        Parameters:

        vals (list): the list of binary values (1/0 or True/False) representing the sample. The input list will not be altered by the operation.

        alpha (float): the significance level for the confidence interval.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(v in {0, 1, True, False} for v in vals):
            raise ValueError("'vals' must contain only binary values (1/0 or True/False).")
        if not isinstance(alpha, float):
            raise TypeError("'alpha' must be a float value.")

        n = len(vals)
        successes = sum(bool(v) for v in vals)
        p_hat = successes / n

        # Standard error
        se = (p_hat * (1 - p_hat) / n) ** 0.5

        # z-critical value (two-tailed)
        z_crit = self._get_critical_value(
            stat_type='z',
            alpha=alpha,
            df=None  # z does not use df
        )

        # Margin of error
        moe = z_crit * se

        lower = p_hat - moe
        upper = p_hat + moe

        return {
            "proportion": round(p_hat, decimals),
            "lower": round(lower, decimals),
            "upper": round(upper, decimals),
            "margin_of_error": round(moe, decimals)
        }

    # ============================================================
    # HYPOTHESIS TESTING - CATEGORICAL
    # ============================================================    
    
    def chi_square(self, vals,  alpha=0.05, decimals=2):
        """Calculate the Chi-squared value for an input list ('vals'). Returns a dictionary containing the chi-squared statistic (key: chi_square), degrees of freedom (key: df), significance decision (key: significance), and the critical value used (key: critical_value).

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        alpha (float): the significance level for the chi-square test.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(vals, list) or not vals:
            raise TypeError("'vals' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All elements in 'vals' must be numeric.")
        if not isinstance(alpha, float):
            raise TypeError("'alpha' must be a float value.")

        n = len(vals)
        df = n - 1

        expected = sum(vals) / n

        chi_sqr = 0
        for v in vals:
            chi_sqr += ((v - expected) ** 2) / expected

        chi_sqr = round(chi_sqr, decimals)

        # Use your internal significance checker
        significance_result = self._check_significance(
            stat_value=chi_sqr,
            df=df,
            test_type='chi2',
            alpha=alpha
        )

        # Extract critical value from the returned string
        # Format: "Likely Significant ... Critical Value: X (DF used: Y)."
        try:
            crit_part = significance_result.split("Critical Value:")[1]
            critical_value = float(crit_part.split("(")[0].strip())
            df_used = crit_part.split("DF used:")[1].replace(").", "").strip()
        except Exception:
            critical_value = None
            df_used = None

        return {
            "chi_square": chi_sqr,
            "df": df,
            "alpha": alpha,
            "critical_value": critical_value,
            "df_used": df_used,
            "significance": significance_result
        }

    def chi_square_independence(self, c_table, alpha=0.05, decimals=2):
        """Perform a Chi-Square Test of Independence on a contingency table represented by a nested list. Returns a dictionary containing the chi-square statistic (key: chi_square), degrees of freedom (key: df), expected frequencies (key: expected), significance decision (key: significance), and the critical value used (key: critical_value).

        Parameters:

        c_table (list): a nested list representing the numeric contingency matrix. This should be the 'c_table' value extracted from the dictionary returned by the contingency_table() function, not the full dictionary itself.

        alpha (float): the significance level for the chi-square test.

        decimals (integer): the number of decimal places float values will be rounded to.
        """
        # Validation
        if not isinstance(c_table, list) or not c_table:
            raise TypeError("'table' must be a non-empty list of lists.")
        if not all(isinstance(row, list) and row for row in c_table):
            raise TypeError("Each row in 'c_table' must be a non-empty list.")
        if len({len(row) for row in c_table}) != 1:
            raise ValueError("All rows in the contingency table must have the same length.")
        if not all(all(isinstance(v, (int, float)) and v >= 0 for v in row) for row in c_table):
            raise ValueError("All values in the contingency table must be non-negative numbers.")
        if not isinstance(alpha, float):
            raise TypeError("'alpha' must be a float value.")

        num_rows = len(c_table)
        num_cols = len(c_table[0])

        # Row totals, column totals, and grand total
        row_totals = [sum(row) for row in c_table]
        col_totals = [sum(c_table[r][c] for r in range(num_rows)) for c in range(num_cols)]
        grand_total = sum(row_totals)

        # Expected frequencies
        expected = [
            [(row_totals[r] * col_totals[c]) / grand_total for c in range(num_cols)]
            for r in range(num_rows)
        ]

        # Chi-square statistic
        chi_sqr = 0
        for r in range(num_rows):
            for c in range(num_cols):
                e = expected[r][c]
                if e > 0:
                    chi_sqr += ((c_table[r][c] - e) ** 2) / e

        chi_sqr = round(chi_sqr, decimals)

        # Degrees of freedom
        df = (num_rows - 1) * (num_cols - 1)

        # Significance check (chi-square is one-tailed)
        significance_result = self._check_significance(
            stat_value=chi_sqr,
            df=df,
            test_type='chi2',
            alpha=alpha
        )

        # Extract critical value and df_used
        try:
            crit_part = significance_result.split("Critical Value:")[1]
            critical_value = float(crit_part.split("(")[0].strip())
            df_used = crit_part.split("DF used:")[1].replace(").", "").strip()
        except Exception:
            critical_value = None
            df_used = None

        # Round expected frequencies
        expected_rounded = [
            [round(e, decimals) for e in row]
            for row in expected
        ]

        return {
            "chi_square": chi_sqr,
            "df": df,
            "alpha": alpha,
            "expected": expected_rounded,
            "critical_value": critical_value,
            "df_used": df_used,
            "significance": significance_result
        }

    def odds_ratio(self, c_table, decimals=2):
        """Calculate the Odds Ratio for a 2x2 contingency table and return a dictionary containing odds_ratio (key: odds_ratio).

        Parameters:

        c_table (list): a nested list representing the numeric contingency matrix. This should be the 'c_table' value extracted from the dictionary returned by the contingency_table() function, not the full dictionary itself.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(c_table, list) or len(c_table) != 2:
            raise ValueError("Odds ratio requires a 2x2 contingency table.")
        if not all(isinstance(row, list) and len(row) == 2 for row in c_table):
            raise ValueError("Odds ratio requires a 2x2 contingency table.")

        a, b = c_table[0]
        c, d = c_table[1]

        if b == 0 or c == 0:
            raise ValueError("Cannot compute odds ratio when b or c is zero.")

        OR = (a * d) / (b * c)
        OR = round(OR, decimals)

        return {"odds_ratio": OR}

    def relative_risk(self, c_table, decimals=2):
        """Calculate the Relative Risk for a 2x2 contingency table and return a dictionary containing relative_risk (key: relative_risk).

        Parameters:

        c_table (list): a nested list representing the numeric contingency matrix. This should be the 'c_table' value extracted from the dictionary returned by the contingency_table() function, not the full dictionary itself.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(c_table, list) or len(c_table) != 2:
            raise ValueError("Relative risk requires a 2x2 contingency table.")
        if not all(isinstance(row, list) and len(row) == 2 for row in c_table):
            raise ValueError("Relative risk requires a 2x2 contingency table.")

        a, b = c_table[0]
        c, d = c_table[1]

        risk1 = a / (a + b)
        risk2 = c / (c + d)

        if risk2 == 0:
            raise ValueError("Cannot compute relative risk when the second risk is zero.")

        RR = risk1 / risk2
        RR = round(RR, decimals)

        return {"relative_risk": RR}

    # ============================================================
    # HYPOTHESIS TESTING - NONPARAMETRIC
    # ============================================================

    def f_test(self, x_var, y_var, alpha=0.05, decimals=2):
        """Perform an F-test to compare the variances of two independent samples ('x_var', 'y_var'). Returns a dictionary containing the F-statistic (key: f_stat), degrees of freedom for both samples (key: df1, key: df2), significance decision (key: significance), and the critical value used (key: critical_value).

        Parameters:

        x_var (list): the list of numeric values representing the first sample. The input list will not be altered by the operation.

        y_var (list): the list of numeric values representing the second sample. The input list will not be altered by the operation.

        alpha (float): the significance level for the F-test.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(x_var, list) or not x_var:
            raise TypeError("'x_var' must be a non-empty list.")
        if not isinstance(y_var, list) or not y_var:
            raise TypeError("'y_var' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in x_var):
            raise ValueError("All elements in 'x_var' must be numeric.")
        if not all(isinstance(v, (int, float)) for v in y_var):
            raise ValueError("All elements in 'y_var' must be numeric.")
        if not isinstance(alpha, float):
            raise TypeError("'alpha' must be a float value.")

        nx = len(x_var)
        ny = len(y_var)

        # Sample means
        xmean = sum(x_var) / nx
        ymean = sum(y_var) / ny

        # Sample variances (unbiased)
        var_x = sum((v - xmean) ** 2 for v in x_var) / (nx - 1)
        var_y = sum((v - ymean) ** 2 for v in y_var) / (ny - 1)

        # Determine which variance is larger (F must be >= 1)
        if var_x >= var_y:
            f_stat = var_x / var_y
            df1 = nx - 1
            df2 = ny - 1
        else:
            f_stat = var_y / var_x
            df1 = ny - 1
            df2 = nx - 1

        f_stat = round(f_stat, decimals)

        # Two-tailed F-test: use 'two_tail' in your significance checker
        significance_result = self._check_significance(
            stat_value=f_stat,
            df=(df1, df2),
            test_type='f',
            alpha=alpha,
            tail_type='two_tail'
        )

        # Extract critical value and df_used
        try:
            crit_part = significance_result.split("Critical Value:")[1]
            critical_value = float(crit_part.split("(")[0].strip())
            df_used = crit_part.split("DF used:")[1].replace(").", "").strip()
        except Exception:
            critical_value = None
            df_used = None

        return {
            "f_stat": f_stat,
            "df1": df1,
            "df2": df2,
            "alpha": alpha,
            "critical_value": critical_value,
            "df_used": df_used,
            "significance": significance_result
        }

    def one_way_anova(self, groups, alpha=0.05, decimals=2):
        """Perform a one-way ANOVA test on k independent samples. Returns a dictionary containing the F-statistic (key: f_stat), degrees of freedom for between-group and within-group variation (key: df_between, key: df_within), significance decision (key: significance), the critical value used (key: critical_value), as well as eta-squared (key: eta_squared) and omega-squared (key: omega_squared).

        Parameters:

        groups (list): a nested list where each element is a list of numeric values representing one independent sample. The input lists will not be altered by the operation.

        alpha (float): the significance level for the ANOVA test.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(groups, list) or len(groups) < 2:
            raise TypeError("'groups' must be a list containing at least two samples.")
        for g in groups:
            if not isinstance(g, list) or not g:
                raise TypeError("Each group must be a non-empty list.")
            if not all(isinstance(v, (int, float)) for v in g):
                raise ValueError("All elements in each group must be numeric.")
        if not isinstance(alpha, float):
            raise TypeError("'alpha' must be a float value.")

        # Number of groups and total sample size
        k = len(groups)
        n_total = sum(len(g) for g in groups)

        # Grand mean
        grand_mean = sum(sum(g) for g in groups) / n_total

        # Between-group sum of squares (SSB)
        ss_between = sum(len(g) * (sum(g) / len(g) - grand_mean) ** 2 for g in groups)

        # Within-group sum of squares (SSW)
        ss_within = 0
        for g in groups:
            mean_g = sum(g) / len(g)
            ss_within += sum((v - mean_g) ** 2 for v in g)

        # Degrees of freedom
        df_between = k - 1
        df_within = n_total - k

        # Mean squares
        ms_between = ss_between / df_between
        ms_within = ss_within / df_within

        # F-statistic
        f_stat = ms_between / ms_within
        f_stat = round(f_stat, decimals)

        # Significance check (F-test, two-tailed)
        significance_result = self._check_significance(
            stat_value=f_stat,
            df=(df_between, df_within),
            test_type='f',
            alpha=alpha,
            tail_type='two_tail'
        )

        # Extract critical value and df_used
        try:
            crit_part = significance_result.split("Critical Value:")[1]
            critical_value = float(crit_part.split("(")[0].strip())
            df_used = crit_part.split("DF used:")[1].replace(").", "").strip()
        except Exception:
            critical_value = None
            df_used = None

        # Total sum of squares
        ss_total = ss_between + ss_within

        # Effect sizes
        eta_squared = ss_between / ss_total if ss_total else 0
        omega_squared = (
            (ss_between - (df_between * ms_within)) /
            (ss_total + ms_within)
        ) if (ss_total + ms_within) else 0

        return {
            "f_stat": f_stat,
            "df_between": df_between,
            "df_within": df_within,
            "alpha": alpha,
            "critical_value": critical_value,
            "df_used": df_used,
            "significance": significance_result,
            "eta_squared": round(eta_squared, decimals),
            "omega_squared": round(omega_squared, decimals)
        }

    def kruskal_wallis(self, groups, alpha=0.05, decimals=2):
        """Perform a Kruskal-Wallis H test on k independent samples. Returns a dictionary containing the H-statistic (key: h_stat), degrees of freedom (key: df), significance decision (key: significance), and the critical value used (key: critical_value).

        Parameters:

        groups (list): a nested list where each element is a list of numeric values representing one independent sample. The input lists will not be altered by the operation.

        alpha (float): the significance level for the Kruskal-Wallis test.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(groups, list) or len(groups) < 2:
            raise TypeError("'groups' must be a list containing at least two samples.")
        for g in groups:
            if not isinstance(g, list) or not g:
                raise TypeError("Each group must be a non-empty list.")
            if not all(isinstance(v, (int, float)) for v in g):
                raise ValueError("All elements in each group must be numeric.")
        if not isinstance(alpha, float):
            raise TypeError("'alpha' must be a float value.")

        # Flatten all values with group labels
        combined = []
        for i, g in enumerate(groups):
            for v in g:
                combined.append((v, i))

        # Sort by value
        combined.sort(key=lambda x: x[0])

        # Assign ranks (with tie handling)
        ranks = [0] * len(combined)
        i = 0
        while i < len(combined):
            j = i
            # Find tie block
            while j < len(combined) and combined[j][0] == combined[i][0]:
                j += 1
            # Average rank for ties
            avg_rank = (i + 1 + j) / 2
            for k in range(i, j):
                ranks[k] = avg_rank
            i = j

        # Sum ranks per group
        rank_sums = [0] * len(groups)
        for (value, group_index), rank in zip(combined, ranks):
            rank_sums[group_index] += rank

        # Total sample size
        n_total = len(combined)

        # Kruskal-Wallis H statistic
        H = 0
        for g, R in zip(groups, rank_sums):
            H += (R ** 2) / len(g)

        H = (12 / (n_total * (n_total + 1))) * H - 3 * (n_total + 1)
        H = round(H, decimals)

        # Degrees of freedom
        df = len(groups) - 1

        # Significance check (chi-square distribution)
        significance_result = self._check_significance(
            stat_value=H,
            df=df,
            test_type='chi2',
            alpha=alpha
        )

        # Extract critical value and df_used
        try:
            crit_part = significance_result.split("Critical Value:")[1]
            critical_value = float(crit_part.split("(")[0].strip())
            df_used = crit_part.split("DF used:")[1].replace(").", "").strip()
        except Exception:
            critical_value = None
            df_used = None

        return {
            "h_stat": H,
            "df": df,
            "alpha": alpha,
            "critical_value": critical_value,
            "df_used": df_used,
            "significance": significance_result
        }

    def mann_whitney(self, x_var, y_var, alpha=0.05, decimals=2):
        """Perform a Mann-Whitney U test on two independent samples ('x_var', 'y_var'). Returns a dictionary containing the U statistic (key: U), the Z statistic (key: z_stat), significance decision (key: significance), and the critical value used (key: critical_value).

        Parameters:

        x_var (list): the list of numeric values representing the first sample. The input list will not be altered by the operation.

        y_var (list): the list of numeric values representing the second sample. The input list will not be altered by the operation.

        alpha (float): the significance level for the Mann-Whitney test.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(x_var, list) or not x_var:
            raise TypeError("'x_var' must be a non-empty list.")
        if not isinstance(y_var, list) or not y_var:
            raise TypeError("'y_var' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in x_var):
            raise ValueError("All elements in 'x_var' must be numeric.")
        if not all(isinstance(v, (int, float)) for v in y_var):
            raise ValueError("All elements in 'y_var' must be numeric.")
        if not isinstance(alpha, float):
            raise TypeError("'alpha' must be a float value.")

        nx = len(x_var)
        ny = len(y_var)

        # Combine values with group labels
        combined = []
        for v in x_var:
            combined.append((v, "x"))
        for v in y_var:
            combined.append((v, "y"))

        # Sort by value
        combined.sort(key=lambda x: x[0])

        # Assign ranks with tie handling
        ranks = [0] * len(combined)
        i = 0
        while i < len(combined):
            j = i
            while j < len(combined) and combined[j][0] == combined[i][0]:
                j += 1
            avg_rank = (i + 1 + j) / 2
            for k in range(i, j):
                ranks[k] = avg_rank
            i = j

        # Sum ranks for each group
        R_x = sum(rank for (val, grp), rank in zip(combined, ranks) if grp == "x")
        R_y = sum(rank for (val, grp), rank in zip(combined, ranks) if grp == "y")

        # Compute U statistics
        U1 = R_x - nx * (nx + 1) / 2
        U2 = R_y - ny * (ny + 1) / 2
        U = min(U1, U2)

        # Normal approximation for Z
        mean_U = nx * ny / 2
        sd_U = math.sqrt(nx * ny * (nx + ny + 1) / 12)

        z_stat = (U - mean_U) / sd_U
        z_stat = round(z_stat, decimals)

        # Use Z-distribution row ('max') in T-table
        significance_result = self._check_significance(
            stat_value=z_stat,
            df=999999,  # triggers Z row
            test_type='t',
            alpha=alpha,
            tail_type='two_tail'
        )

        # Extract critical value and df_used
        try:
            crit_part = significance_result.split("Critical Value:")[1]
            critical_value = float(crit_part.split("(")[0].strip())
            df_used = crit_part.split("DF used:")[1].replace(").", "").strip()
        except Exception:
            critical_value = None
            df_used = None

        return {
            "U": round(U, decimals),
            "z_stat": z_stat,
            "nx": nx,
            "ny": ny,
            "alpha": alpha,
            "critical_value": critical_value,
            "df_used": df_used,
            "significance": significance_result
        }

    def wilcoxon_signed_rank(self, x_var, y_var, alpha=0.05, decimals=2):
        """Perform a Wilcoxon Signed-Rank Test on two paired samples ('x_var', 'y_var'). Returns a dictionary containing the W statistic (key: W), the Z statistic (key: z_stat), significance decision (key: significance), and the critical value used (key: critical_value).

        Parameters:

        x_var (list): the list of numeric values representing the first paired sample. The input list will not be altered by the operation.

        y_var (list): the list of numeric values representing the second paired sample. The input list will not be altered by the operation.

        alpha (float): the significance level for the Wilcoxon test.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(x_var, list) or not x_var:
            raise TypeError("'x_var' must be a non-empty list.")
        if not isinstance(y_var, list) or not y_var:
            raise TypeError("'y_var' must be a non-empty list.")
        if len(x_var) != len(y_var):
            raise ValueError("'x_var' and 'y_var' must be the same length.")
        if not all(isinstance(v, (int, float)) for v in x_var):
            raise ValueError("All elements in 'x_var' must be numeric.")
        if not all(isinstance(v, (int, float)) for v in y_var):
            raise ValueError("All elements in 'y_var' must be numeric.")
        if not isinstance(alpha, float):
            raise TypeError("'alpha' must be a float value.")

        # Compute paired differences
        diffs = [y - x for x, y in zip(x_var, y_var)]

        # Remove zero differences (they contribute no rank)
        nonzero = [(abs(d), d) for d in diffs if d != 0]

        if not nonzero:
            raise ValueError("All paired differences are zero; Wilcoxon test cannot be performed.")

        # Sort by absolute difference
        nonzero.sort(key=lambda x: x[0])

        # Assign ranks with tie handling
        ranks = [0] * len(nonzero)
        i = 0
        while i < len(nonzero):
            j = i
            while j < len(nonzero) and nonzero[j][0] == nonzero[i][0]:
                j += 1
            avg_rank = (i + 1 + j) / 2
            for k in range(i, j):
                ranks[k] = avg_rank
            i = j

        # Compute W+ and W- (sum of signed ranks)
        W_pos = sum(rank for rank, (_, d) in zip(ranks, nonzero) if d > 0)
        W_neg = sum(rank for rank, (_, d) in zip(ranks, nonzero) if d < 0)

        # Test statistic W is the smaller of W+ and W-
        W = min(W_pos, W_neg)
        W = round(W, decimals)

        # Normal approximation for Z
        n = len(nonzero)
        mean_W = n * (n + 1) / 4
        sd_W = math.sqrt(n * (n + 1) * (2*n + 1) / 24)

        z_stat = (W - mean_W) / sd_W
        z_stat = round(z_stat, decimals)

        # Use Z-distribution row ('max') in T-table
        significance_result = self._check_significance(
            stat_value=z_stat,
            df=999999,  # triggers Z row
            test_type='t',
            alpha=alpha,
            tail_type='two_tail'
        )

        # Extract critical value and df_used
        try:
            crit_part = significance_result.split("Critical Value:")[1]
            critical_value = float(crit_part.split("(")[0].strip())
            df_used = crit_part.split("DF used:")[1].replace(").", "").strip()
        except Exception:
            critical_value = None
            df_used = None

        return {
            "W": W,
            "z_stat": z_stat,
            "n": n,
            "alpha": alpha,
            "critical_value": critical_value,
            "df_used": df_used,
            "significance": significance_result
        }

    # ============================================================
    # EFFECT SIZES
    # ============================================================

    def effect_size_d_g(self, group1, group2, decimals=2):
        """Calculate and return Cohen's d and Hedges' g for two independent groups ('group1' and 'group2').  Returns a dictionary with the following keys: cohens_d and hedges_g.

        Parameters:

        group1 (list): the list of numeric values representing the first group. The input list will not be altered by the operation.

        group2 (list): the list of numeric values representing the second group. The input list will not be altered by the operation.

        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validation
        if not isinstance(group1, list) or not group1:
            raise TypeError("'group1' must be a non-empty list.")
        if not isinstance(group2, list) or not group2:
            raise TypeError("'group2' must be a non-empty list.")
        if not all(isinstance(v, (int, float)) for v in group1):
            raise ValueError("'group1' must contain only numeric values.")
        if not all(isinstance(v, (int, float)) for v in group2):
            raise ValueError("'group2' must contain only numeric values.")

        n1 = len(group1)
        n2 = len(group2)

        mean1 = sum(group1) / n1
        mean2 = sum(group2) / n2

        # Sample variances
        var1 = sum((x - mean1)**2 for x in group1) / (n1 - 1)
        var2 = sum((x - mean2)**2 for x in group2) / (n2 - 1)

        # Pooled standard deviation
        sp = (( (n1 - 1)*var1 + (n2 - 1)*var2 ) / (n1 + n2 - 2)) ** 0.5

        # Cohen's d
        d = (mean1 - mean2) / sp if sp else 0

        # Hedges' correction factor
        J = 1 - (3 / (4*(n1 + n2) - 9))

        # Hedges' g
        g = J * d

        return {
            "cohens_d": round(d, decimals),
            "hedges_g": round(g, decimals)
        }
  
class _Conversion():
    """A set of functions for performing unit conversions on lists of numerical values."""
    # ============================================================
    # PRIVATE METHODS AND CLASS ATTRIBUTES AND FUNCTION DISPLAY
    # ============================================================
    def __init__(self):
        self._total_funcs = len([f for f in dir(__class__) if not f.startswith('__')])
    
    def __repr__(self):
        """If the toolbox is printed, display a message."""
        return ".conversion: a set of functions for performing unit conversions on lists of numerical values."
    
    def display_toolbox_functions(self):
        """Display a list of all available functions within this toolbox."""
        print(f"Number of {__class__.__name__[1:]} functions: {len([f for f in dir(__class__) if not f.startswith('__')])}")
        for f in [f for f in dir(__class__) if not f.startswith("__")]:
            print(f)
    
    # ============================================================
    # LENGTH, AREA, VOLUME
    # ============================================================
    
    def convert_distance(self, vals, from_unit, to_unit, decimals=2):
        """Convert the values of the input list ('vals') from their current unit of distance ('from_unit') to another unit of distance ('to_unit') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        from_unit (string): the unit the list of values is being converted from. Possible values include "mm" for millimeters, "cm" for centimeters, "m" for meters, "km" for kilometers, "in" for inches, "ft" for feet, "yd" for yards, and "mi" for miles.

        to_unit (string): the unit the list of values is being converted to. Possible values include "mm" for millimeters, "cm" for centimeters, "m" for meters, "km" for kilometers, "in" for inches, "ft" for feet, "yd" for yards, and "mi" for miles.

        decimals (integer): the number of decimal places float values will be rounded to."""

        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All values must be numeric.")

        units = {"mm": 0.001,
                 "cm": 0.01,
                 "m": 1.0,
                 "km": 1000.0,
                 "in": 0.0254,
                 "ft": 0.3048,
                 "yd": 0.9144,
                 "mi": 1609.344}

        if from_unit not in units or to_unit not in units:
            raise ValueError("Unsupported unit.")

        # Conversion via a reference scale
        factor_from = units[from_unit]
        factor_to = units[to_unit]

        result = [(v * factor_from) / factor_to for v in vals]
        return [round(r, decimals) for r in result]

    def convert_area(self, vals, from_unit, to_unit, decimals=2):
        """Convert the values of the input list ('vals') from their current unit of area ('from_unit') to another unit of area ('to_unit') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        from_unit (string): the unit the list of values is being converted from. Possible values include "cm" for centimeters squared, "m" for meters squared, "km" for kilometers squared, "ha" for hectares, "ft" for feet squared, "mi" for miles squared, and "ac" for acres.

        to_unit (string): the unit the list of values is being converted to. Possible values include "cm" for centimeters squared, "m" for meters squared, "km" for kilometers squared, "ha" for hectares, "ft" for feet squared, "mi" for miles squared, and "ac" for acres.

        decimals (integer): the number of decimal places float values will be rounded to."""

        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All values must be numeric.")

        units = {
            "cm": 0.0001,          # 1 cm² = 0.0001 reference units
            "m": 1.0,              # 1 m² = 1 reference unit
            "km": 1_000_000.0,     # 1 km² = 1,000,000 reference units
            "ha": 10_000.0,        # 1 hectare = 10,000 reference units
            "ft": 0.09290304,      # 1 ft² = 0.09290304 reference units
            "mi": 2_589_988.11,    # 1 mi² = 2,589,988.11 reference units
            "ac": 4046.8564224     # 1 acre = 4046.8564224 reference units
        }

        if from_unit not in units or to_unit not in units:
            raise ValueError("Unsupported unit.")

        factor_from = units[from_unit]
        factor_to = units[to_unit]

        result = [(v * factor_from) / factor_to for v in vals]
        return [round(r, decimals) for r in result]

    def convert_volume(self, vals, from_unit, to_unit, decimals=2):
        """Convert the values of the input list ('vals') from their current unit of volume ('from_unit') to another unit of volume ('to_unit') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        from_unit (string): the unit the list of values is being converted from. Possible values include "ml" for millilitres, "l" for litres, "mm" for cubic millimeters, "cm" for cubic centimeters, "m" for cubic meters, "in" for cubic inches, "ft" for cubic feet, "oz" for fluid ounces, and "gal" for gallons.

        to_unit (string): the unit the list of values is being converted to. Possible values include "ml" for millilitres, "l" for litres, "mm" for cubic millimeters, "cm" for cubic centimeters, "m" for cubic meters, "in" for cubic inches, "ft" for cubic feet, "oz" for fluid ounces, and "gal" for gallons.

        decimals (integer): the number of decimal places float values will be rounded to."""

        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All values must be numeric.")

        units = {
            "ml": 0.001,             # 1 millilitre = 0.001 reference units
            "l": 1.0,                # 1 litre = 1 reference unit
            "mm": 1e-9,              # 1 cubic millimetre = 1e-9 reference units
            "cm": 1e-6,              # 1 cubic centimetre = 1e-6 reference units
            "m": 1000.0,             # 1 cubic metre = 1000 reference units
            "in": 0.016387064,       # 1 cubic inch = 0.016387064 reference units
            "ft": 28.316846592,      # 1 cubic foot = 28.316846592 reference units
            "oz": 0.0295735295625,   # 1 US fluid ounce = 0.0295735295625 reference units
            "gal": 3.785411784       # 1 US gallon = 3.785411784 reference units
        }

        if from_unit not in units or to_unit not in units:
            raise ValueError("Unsupported unit.")

        factor_from = units[from_unit]
        factor_to = units[to_unit]

        result = [(v * factor_from) / factor_to for v in vals]
        return [round(r, decimals) for r in result]

    # ============================================================
    # MASS AND WEIGHT
    # ============================================================

    def convert_weight(self, vals, from_unit, to_unit, decimals=2):
        """Convert the values of the input list ('vals') from their current unit of weight ('from_unit') to another unit of weight ('to_unit') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        from_unit (string): the unit the list of values is being converted from. Possible values include "mg" for milligrams, "g" for grams, "kg" for kilograms, "t" for metric tonnes, "oz" for ounces, "lb" for pounds, "st" for imperial short tons, and "lt" for imperial long tons.

        to_unit (string): the unit the list of values is being converted to. Possible values include "mg" for milligrams, "g" for grams, "kg" for kilograms, "t" for metric tonnes, "oz" for ounces, "lb" for pounds, "st" for imperial short tons, and "lt" for imperial long tons.

        decimals (integer): the number of decimal places float values will be rounded to."""

        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All values must be numeric.")

        units = {
            "mg": 0.001,             # 1 milligram = 0.001 reference units
            "g": 1.0,                # 1 gram = 1 reference unit
            "kg": 1000.0,            # 1 kilogram = 1000 reference units
            "t": 1_000_000.0,        # 1 metric tonne = 1,000,000 reference units
            "oz": 28.349523125,      # 1 ounce = 28.349523125 reference units
            "lb": 453.59237,         # 1 pound = 453.59237 reference units
            "st": 907_184.74,        # 1 short ton = 907,184.74 reference units
            "lt": 1_016_046.9088     # 1 long ton = 1,016,046.9088 reference units
        }

        if from_unit not in units or to_unit not in units:
            raise ValueError("Unsupported unit.")

        factor_from = units[from_unit]
        factor_to = units[to_unit]

        result = [(v * factor_from) / factor_to for v in vals]
        return [round(r, decimals) for r in result]

    def convert_density(self, vals, from_unit, to_unit, decimals=2):
        """Convert the values of the input list ('vals') from their current unit of density ('from_unit') to another unit of density ('to_unit') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        from_unit (string): the unit the list of values is being converted from. Possible values include "kg/m3" for kilograms per cubic meter, "g/cm3" for grams per cubic centimeter, and "lb/ft3" for pounds per cubic foot.

        to_unit (string): the unit the list of values is being converted to. Possible values include "kg/m3" for kilograms per cubic meter, "g/cm3" for grams per cubic centimeter, and "lb/ft3" for pounds per cubic foot.

        decimals (integer): the number of decimal places float values will be rounded to."""

        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All values must be numeric.")

        units = {
            "kg/m3": 1.0,             # 1 kilogram per cubic meter = 1 reference unit
            "g/cm3": 1000.0,          # 1 gram per cubic centimeter = 1000 reference units
            "lb/ft3": 16.0185         # 1 pound per cubic foot = 16.0185 reference units
        }

        if from_unit not in units or to_unit not in units:
            raise ValueError("Unsupported unit.")

        factor_from = units[from_unit]
        factor_to = units[to_unit]

        result = [(v * factor_from) / factor_to for v in vals]
        return [round(r, decimals) for r in result]

    # ============================================================
    # TIME AND VELOCITY
    # ============================================================

    def convert_time(self, vals, from_unit, to_unit, decimals=2):
        """Convert the values of the input list ('vals') from their current unit of time ('from_unit') to another unit of time ('to_unit') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        from_unit (string): the unit the list of values is being converted from. Possible values include "ms" for milliseconds, "s" for seconds, "m" for minutes, "hr" for hours, "d" for days, and "y" for years.

        to_unit (string): the unit the list of values is being converted to. Possible values include "ms" for milliseconds, "s" for seconds, "m" for minutes, "hr" for hours, "d" for days, and "y" for years.

        decimals (integer): the number of decimal places float values will be rounded to."""

        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All values must be numeric.")

        units = {
            "ms": 0.001,             # 1 millisecond = 0.001 reference units
            "s": 1.0,                # 1 second = 1 reference unit
            "m": 60.0,               # 1 minute = 60 reference units
            "hr": 3600.0,            # 1 hour = 3600 reference units
            "d": 86400.0,            # 1 day = 86,400 reference units
            "y": 31_556_952.0        # 1 year (365.2425 days) = 31,556,952 reference units
        }

        if from_unit not in units or to_unit not in units:
            raise ValueError("Unsupported unit.")

        factor_from = units[from_unit]
        factor_to = units[to_unit]

        result = [(v * factor_from) / factor_to for v in vals]
        return [round(r, decimals) for r in result]

    def convert_velocity(self, vals, from_unit, to_unit, decimals=2):
        """Convert the values of the input list ('vals') from their current unit of velocity ('from_unit') to another unit of velocity ('to_unit') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        from_unit (string): the unit the list of values is being converted from. Possible values include "m/s" for meters per second, "km/h" for kilometers per hour, "mph" for miles per hour, "ft/s" for feet per second, and "kn" for knots.

        to_unit (string): the unit the list of values is being converted to. Possible values include "m/s" for meters per second, "km/h" for kilometers per hour, "mph" for miles per hour, "ft/s" for feet per second, and "kn" for knots.

        decimals (integer): the number of decimal places float values will be rounded to."""

        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All values must be numeric.")

        units = {
            "m/s": 1.0,             # 1 meter per second = 1 reference unit
            "km/h": 0.277778,       # 1 kilometer per hour = 0.277778 reference units
            "mph": 0.44704,         # 1 mile per hour = 0.44704 reference units
            "ft/s": 0.3048,         # 1 foot per second = 0.3048 reference units
            "kn": 0.514444          # 1 knot = 0.514444 reference units
        }

        if from_unit not in units or to_unit not in units:
            raise ValueError("Unsupported unit.")

        factor_from = units[from_unit]
        factor_to = units[to_unit]

        result = [(v * factor_from) / factor_to for v in vals]
        return [round(r, decimals) for r in result]

    def convert_flowrate(self, vals, from_unit, to_unit, decimals=2):
        """Convert the values of the input list ('vals') from their current unit of flow rate ('from_unit') to another unit of flow rate ('to_unit') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        from_unit (string): the unit the list of values is being converted from. Possible values include "m3/s" for cubic meters per second, "L/s" for liters per second, "L/min" for liters per minute, "gal/min" for US gallons per minute, and "ft3/s" for cubic feet per second.

        to_unit (string): the unit the list of values is being converted to. Possible values include "m3/s" for cubic meters per second, "L/s" for liters per second, "L/min" for liters per minute, "gal/min" for US gallons per minute, and "ft3/s" for cubic feet per second.

        decimals (integer): the number of decimal places float values will be rounded to."""

        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All values must be numeric.")

        units = {
            "m3/s": 1.0,              # 1 cubic meter per second = 1 reference unit
            "L/s": 0.001,             # 1 liter per second = 0.001 reference units
            "L/min": 0.001 / 60.0,    # 1 liter per minute = 0.001/60 reference units
            "gal/min": 0.00378541 / 60.0,  # 1 US gallon per minute = 0.00378541/60 reference units
            "ft3/s": 0.0283168        # 1 cubic foot per second = 0.0283168 reference units
        }

        if from_unit not in units or to_unit not in units:
            raise ValueError("Unsupported unit.")

        factor_from = units[from_unit]
        factor_to = units[to_unit]

        result = [(v * factor_from) / factor_to for v in vals]
        return [round(r, decimals) for r in result]

    # ============================================================
    # PRESSURE, ENERGY, POWER, FORCE
    # ============================================================
    
    def convert_pressure(self, vals, from_unit, to_unit, decimals=2):
        """Convert the values of the input list ('vals') from their current unit of pressure ('from_unit') to another unit of pressure ('to_unit') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        from_unit (string): the unit the list of values is being converted from. Possible values include "Pa" for pascals, "kPa" for kilopascals, "bar" for bars, "atm" for atmospheres, "psi" for pounds per square inch, and "mmHg" for millimeters of mercury.

        to_unit (string): the unit the list of values is being converted to. Possible values include "Pa" for pascals, "kPa" for kilopascals, "bar" for bars, "atm" for atmospheres, "psi" for pounds per square inch, and "mmHg" for millimeters of mercury.

        decimals (integer): the number of decimal places float values will be rounded to."""

        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All values must be numeric.")

        units = {
            "Pa": 1.0,               # 1 pascal = 1 reference unit
            "kPa": 1000.0,           # 1 kilopascal = 1000 reference units
            "bar": 100000.0,         # 1 bar = 100,000 reference units
            "atm": 101325.0,         # 1 atmosphere = 101,325 reference units
            "psi": 6894.757,         # 1 psi = 6,894.757 reference units
            "mmHg": 133.322          # 1 mmHg = 133.322 reference units
        }

        if from_unit not in units or to_unit not in units:
            raise ValueError("Unsupported unit.")

        factor_from = units[from_unit]
        factor_to = units[to_unit]

        result = [(v * factor_from) / factor_to for v in vals]
        return [round(r, decimals) for r in result]

    def convert_energy(self, vals, from_unit, to_unit, decimals=2):
        """Convert the values of the input list ('vals') from their current unit of energy ('from_unit') to another unit of energy ('to_unit') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        from_unit (string): the unit the list of values is being converted from. Possible values include "J" for joules, "kJ" for kilojoules, "cal" for calories, "kcal" for kilocalories, "Wh" for watt-hours, "kWh" for kilowatt-hours, and "BTU" for British thermal units.

        to_unit (string): the unit the list of values is being converted to. Possible values include "J" for joules, "kJ" for kilojoules, "cal" for calories, "kcal" for kilocalories, "Wh" for watt-hours, "kWh" for kilowatt-hours, and "BTU" for British thermal units.

        decimals (integer): the number of decimal places float values will be rounded to."""

        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All values must be numeric.")

        units = {
            "J": 1.0,                # 1 joule = 1 reference unit
            "kJ": 1000.0,            # 1 kilojoule = 1000 reference units
            "cal": 4.184,            # 1 calorie (thermochemical) = 4.184 reference units
            "kcal": 4184.0,          # 1 kilocalorie = 4184 reference units
            "Wh": 3600.0,            # 1 watt-hour = 3600 reference units
            "kWh": 3_600_000.0,      # 1 kilowatt-hour = 3,600,000 reference units
            "BTU": 1055.06           # 1 BTU = 1055.06 reference units
        }

        if from_unit not in units or to_unit not in units:
            raise ValueError("Unsupported unit.")

        factor_from = units[from_unit]
        factor_to = units[to_unit]

        result = [(v * factor_from) / factor_to for v in vals]
        return [round(r, decimals) for r in result]

    def convert_power(self, vals, from_unit, to_unit, decimals=2):
        """Convert the values of the input list ('vals') from their current unit of power ('from_unit') to another unit of power ('to_unit') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        from_unit (string): the unit the list of values is being converted from. Possible values include "W" for watts, "kW" for kilowatts, "MW" for megawatts, and "hp" for horsepower.

        to_unit (string): the unit the list of values is being converted to. Possible values include "W" for watts, "kW" for kilowatts, "MW" for megawatts, and "hp" for horsepower.

        decimals (integer): the number of decimal places float values will be rounded to."""

        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All values must be numeric.")

        units = {
            "W": 1.0,                # 1 watt = 1 reference unit
            "kW": 1000.0,            # 1 kilowatt = 1000 reference units
            "MW": 1_000_000.0,       # 1 megawatt = 1,000,000 reference units
            "hp": 745.7              # 1 horsepower (mechanical) = 745.7 reference units
        }

        if from_unit not in units or to_unit not in units:
            raise ValueError("Unsupported unit.")

        factor_from = units[from_unit]
        factor_to = units[to_unit]

        result = [(v * factor_from) / factor_to for v in vals]
        return [round(r, decimals) for r in result]

    def convert_force(self, vals, from_unit, to_unit, decimals=2):
        """Convert the values of the input list ('vals') from their current unit of force ('from_unit') to another unit of force ('to_unit') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        from_unit (string): the unit the list of values is being converted from. Possible values include "N" for newtons, "kN" for kilonewtons, "dyn" for dynes, and "lbf" for pounds-force.

        to_unit (string): the unit the list of values is being converted to. Possible values include "N" for newtons, "kN" for kilonewtons, "dyn" for dynes, and "lbf" for pounds-force.

        decimals (integer): the number of decimal places float values will be rounded to."""

        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All values must be numeric.")

        units = {
            "N": 1.0,                 # 1 newton = 1 reference unit
            "kN": 1000.0,             # 1 kilonewton = 1000 reference units
            "dyn": 1e-5,              # 1 dyne = 1e-5 reference units
            "lbf": 4.44822            # 1 pound-force = 4.44822 reference units
        }

        if from_unit not in units or to_unit not in units:
            raise ValueError("Unsupported unit.")

        factor_from = units[from_unit]
        factor_to = units[to_unit]

        result = [(v * factor_from) / factor_to for v in vals]
        return [round(r, decimals) for r in result]

    # ============================================================
    # TEMPERATURE AND CONCENTRATION
    # ============================================================

    def convert_temperature(self, vals, from_unit, to_unit, decimals=2):
        """Convert the values of the input list ('vals') from their current unit of temperature ('from_unit') to another unit of temperature ('to_unit') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        from_unit (string): the unit the list of values is being converted from. Possible values include "k" for kelvin, "c" for celsius, and "f" for fahrenheit.

        to_unit (string): the unit the list of values is being converted to. Possible values include "k" for kelvin, "c" for celsius, and "f" for fahrenheit.

        decimals (integer): the number of decimal places float values will be rounded to."""

        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All values must be numeric.")

        valid_units = {"k", "c", "f"}
        if from_unit not in valid_units or to_unit not in valid_units:
            raise ValueError("Unsupported unit.")

        # Conversion functions to and from a reference scale (Celsius)
        def to_celsius(value, unit):
            if unit == "c":
                return value
            elif unit == "k":
                return value - 273.15
            elif unit == "f":
                return (value - 32) * 5.0 / 9.0

        def from_celsius(value, unit):
            if unit == "c":
                return value
            elif unit == "k":
                return value + 273.15
            elif unit == "f":
                return (value * 9.0 / 5.0) + 32

        # Convert each value
        result = []
        for v in vals:
            celsius_val = to_celsius(v, from_unit)
            converted = from_celsius(celsius_val, to_unit)
            result.append(round(converted, decimals))

        return result
    
    def convert_concentration(self, vals, from_unit, to_unit, decimals=2):
        """Convert the values of the input list ('vals') from their current unit of concentration ('from_unit') to another unit of concentration ('to_unit') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        from_unit (string): the unit the list of values is being converted from. Possible values include "mg/L" for milligrams per liter, "g/L" for grams per liter, "ppm" for parts per million, and "%" for percent by mass/volume (assuming dilute aqueous solutions).

        to_unit (string): the unit the list of values is being converted to. Possible values include "mg/L" for milligrams per liter, "g/L" for grams per liter, "ppm" for parts per million, and "%" for percent by mass/volume (assuming dilute aqueous solutions).

        decimals (integer): the number of decimal places float values will be rounded to."""

        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not all(isinstance(v, (int, float)) for v in vals):
            raise ValueError("All values must be numeric.")

        units = {
            "mg/L": 0.001,       # 1 mg/L = 0.001 reference units
            "g/L": 1.0,          # 1 g/L = 1 reference unit
            "ppm": 0.001,        # 1 ppm ≈ 1 mg/L = 0.001 reference units (for water-like density)
            "%": 10.0            # 1% = 10 g/L (assuming 1% = 10 g/L in dilute aqueous solutions)
        }

        if from_unit not in units or to_unit not in units:
            raise ValueError("Unsupported unit.")

        factor_from = units[from_unit]
        factor_to = units[to_unit]

        result = [(v * factor_from) / factor_to for v in vals]
        return [round(r, decimals) for r in result]
       
class _Points():
    """A set of functions for performing operations on PointDatasets"""
    # ============================================================
    # PRIVATE METHODS AND CLASS ATTRIBUTES AND FUNCTION DISPLAY
    # ============================================================
    def __init__(self):
        self._total_funcs = len([f for f in dir(__class__) if not f.startswith('__')])
    
    def __repr__(self):
        """If the toolbox is printed, display a message."""
        return ".points: a set of functions for performing operations on PointDatasets."

    def display_toolbox_functions(self):
        """Display a list of all available functions within this toolbox."""
        print(f"Number of {__class__.__name__[1:]} functions: {len([f for f in dir(__class__) if not f.startswith('__')])}")
        for f in [f for f in dir(__class__) if not f.startswith("__")]:
            print(f)

    # ============================================================
    # TABLE CONVERSION
    # ============================================================

    def table_to_points(self, table, id_col, x_col, y_col):
        """Convert a Table object ('table') to a PointDataset containing Point objects. A column must be chosen for the unique identifier ('id_col'), x coordinate ('x_col'), and y coordinate ('y_col'). All other columns will be added as point attributes which are referenced by their column header.

        Parameters:
        
        table (object): the Table object to convert to a PointDataset.
        
        id_col (string): the header of the column that contains the unique identifier for each point.
        
        x_col (string): the header of the column that contains the x coordinate for each point.
        
        y_col (string): the header of the column that contains the y coordinate for each point."""
        if not isinstance(table,_Table):
            raise ValueError("'table' must be a Table object.")
        headers = table.get_headers()
        nc = table.num_columns()

        # Validate required columns exist
        try:
            id_idx = headers.index(id_col)
            x_idx = headers.index(x_col)
            y_idx = headers.index(y_col)
        except ValueError:
            raise ValueError("Specified id/x/y column not found in headers.")
    
        ids = [table.get_row(r)[id_idx] for r in range(table.num_rows())]
        if len(set(ids)) != len(ids):
            raise ValueError(f"Duplicate IDs found in column '{id_col}'. IDs must be unique.")


        points = []
        for r in range(table.num_rows()):
            row = table.get_row(r)
            id_val = row[id_idx]
            x_val = row[x_idx]
            y_val = row[y_idx]

            # Collect all other attributes
            attrs = [row[c] for c in range(nc) if c not in (id_idx, x_idx, y_idx)]
            points.append(_Point(id_val, x_val, y_val, attrs))

        return _PointDataset(points, headers, id_col=id_col, x_col=x_col, y_col=y_col)

    def points_to_table(self, points):
        """Convert a PointDataset ('points') into a Table object.
        
        Parameters:
        
        points (object): the PointDataset to convert to a Table object."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")

        output = TableTools().new_table()
        headers = points.headers

        for p in points:
            # Build row: id, x, y, followed by attrs
            row = [p.id, p.x, p.y] + list(p.attrs)
            output.append_row(row)

        output.set_headers(headers)
        output.detect_dtypes()
        return output

    # ============================================================
    # POINT CREATION AND ACCESS
    # ============================================================

    def get_point(self, points, id):
        """Return a single Point object with the specified id ('id') from a PointDataset ('points').
        
        Parameters:
        
        point_dataset (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        id (integer, float, string): the unique identifier of the point to be returned."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'point_dataset' must be a PointDataset object.")
        ids = [p.id for p in points]
        if id not in ids:
            raise ValueError(f"Point with id '{id}' not found in dataset.")


        for p in points:
            if p.id == id:
                # Return a shallow copy of the Point to avoid mutating the original
                return _Point(p.id, p.x, p.y, list(p.attrs))

    def new_point(self, id, x, y, attr_name=None, attr=None):
        """Return a new Point object with the specified unique identifier, x coordinate, y coordinate, and attributes.
        
        Parameters:
        
        id (integer, float, string): the unique identifier of the point.
        
        x (integer, float): the x coordinate of the point.
        
        y (integer, float): the y coordinate of the point.
        
        attr_name (string, list): a string or list of strings representing the names of point attributes.
        
        attr (integer, float, string, list): the value or list of values representing the values of point attributes."""
        if not isinstance(id, (int, float, str)):
            raise TypeError("'id' must be an integer, float, or string.")
        
        # Validate coordinates
        if not isinstance(x, (int, float)):
            raise TypeError("'x' must be an integer or float.")
        if not isinstance(y, (int, float)):
            raise TypeError("'y' must be an integer or float.")

        if attr_name is not None:
            if isinstance(attr_name, str):
                attr_name = [attr_name]
        if attr is not None:
            if not isinstance(attr, list):
                attr = [attr]

        attrs = []
        if attr_name is not None and attr is not None:
            if len(attr_name) != len(attr):
                raise ValueError("Length of 'attr_name' and 'attr' must match.")
            # Store values only; headers are tracked at the PointDataset level
            attrs = list(attr)

        return _Point(id, x, y, attrs)

    def add_point(self, points, point):
        """Add a new Point object ('point') to a PointDataset ('points').
        
        Parameters:
        
        points (object): the PointDataset to add the point to.
        
        point (object): the Point object to be added."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if not isinstance(point, _Point):
            raise TypeError("'point' must be a Point object.")

        # Validate headers length matches point.attrs length
        expected_attr_count = len(points.headers) - 3  # subtract id, x, y
        if len(point.attrs) != expected_attr_count:
            raise ValueError("Point attributes length does not match dataset headers.")
        
        existing_ids = {p.id for p in points.points}
        if point.id in existing_ids:
            raise ValueError(f"Point ID '{point.id}' already exists in dataset. IDs must be unique.")

        # Add point
        points.points.append(point)

        # Update extent
        points._update_extent()
  
    def sample_points(self, points, num_points):
        """Return a random sample of points from a PointDataset ('points') and return a new PointDataset.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        num_points (integer): the number of points to randomly sample from the dataset."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if not isinstance(num_points, int) or num_points <= 0:
            raise ValueError("'num_points' must be a positive integer.")
        if num_points > len(points.points):
            raise ValueError("'num_points' cannot be greater than the number of points in the dataset.")

        sampled_points = random.sample(points.points, num_points)

        # Copy points to avoid mutating originals
        new_points = [_Point(p.id, p.x, p.y, p.attrs.copy()) for p in sampled_points]

        return _PointDataset(new_points, points.headers.copy(),id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    # ============================================================
    # POINT REMOVAL
    # ============================================================

    def remove_point(self, points, id):
        """Remove a Point object with the specified id ('id') from a PointDataset ('points').
        
        Parameters:
        
        points (object): the PointDataset to remove the point from.
        
        id (integer, float, string): the unique identifier of the point to be removed."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        
        # Validate that the id exists in the dataset
        ids = [p.id for p in points]
        if id not in ids:
            raise ValueError(f"Point with id '{id}' not found in dataset.")
        
        # Remove the point
        points.points = [p for p in points if p.id != id]
        
        # Update extent
        points._update_extent()

    def remove_overlapping_points(self, points, dist_thresh):
        """Remove overlapping points from a PointDataset ('points') and return a new PointDataset.
        
        Parameters:
        
        points (object): the dataset to process. The input dataset will not be altered by the operation.
        
        dist_thresh (int or float): the distance threshold to remove points by. Any point within the distance threshold of another point will be removed."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if not isinstance(dist_thresh, (int, float)):
            raise ValueError("'dist_thresh' must be a number.")

        output_points = []
        for p in points.points:
            remove = False
            for p2 in points.points:
                if p is not p2:  # avoid self-comparison
                    dx = p.x - p2.x
                    dy = p.y - p2.y
                    dist = math.sqrt(dx * dx + dy * dy)
                    if dist <= dist_thresh:
                        remove = True
                        break  # no need to check further
            if not remove:
                output_points.append(_Point(p.id, p.x, p.y, p.attrs.copy()))

        return _PointDataset(output_points, points.headers.copy(), id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)
    
    # ============================================================
    # ATTRIBUTE MANAGEMENT
    # ============================================================

    def add_attribute(self, points, attr_name, attr):
        """Add a new attribute to each Point in a PointDataset ('points') and return a new PointDataset.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        attr_name (string): the name of the attribute to add.
        
        attr (list): the list of values representing the values of the point attribute to be added."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if not isinstance(attr_name, str):
            raise TypeError("'attr_name' must be a string.")
        if attr_name in ("id", "x", "y"):
            raise ValueError(f"Cannot add core attribute '{attr_name}'.")
        if not isinstance(attr, list):
            raise TypeError("'attr' must be a list.")
        if len(attr) != len(points):
            raise ValueError("Length of 'attr' must match number of points in dataset.")

        # Copy points to avoid mutating original dataset
        new_points = []
        for i, p in enumerate(points):
            new_attrs = list(p.attrs) + [attr[i]]
            new_points.append(_Point(p.id, p.x, p.y, new_attrs))

        # Update headers: add new attribute name
        new_headers = list(points.headers) + [attr_name]

        return _PointDataset(new_points, new_headers, id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)
    
    def remove_attribute(self, points, attr_name):
        """Remove an attribute from each Point in a PointDataset ('points') and return a new PointDataset.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        attr_name (string): the name of the attribute to remove."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if not isinstance(attr_name, str):
            raise TypeError("'attr_name' must be a string.")
        if attr_name not in points.headers:
            raise ValueError(f"Attribute '{attr_name}' not found in dataset headers.")
        if attr_name in ("id", "x", "y"):
            raise ValueError(f"Cannot remove core attribute '{attr_name}'.")

        # Find index of attribute to remove
        attr_index = points.headers.index(attr_name) - 3  # offset for id, x, y

        # Copy points with attribute removed
        new_points = []
        for p in points:
            new_attrs = list(p.attrs)
            if 0 <= attr_index < len(new_attrs):
                new_attrs.pop(attr_index)
            new_points.append(_Point(p.id, p.x, p.y, new_attrs))

        # Update headers
        new_headers = [h for h in points.headers if h != attr_name]

        return _PointDataset(new_points, new_headers, id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    def rename_attribute(self, points, old_name, new_name):
        """Rename an attribute in a PointDataset ('points') and return a new PointDataset.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        old_name (string): the current name of the attribute to rename.
        
        new_name (string): the new name for the attribute."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if not isinstance(old_name, str) or not isinstance(new_name, str):
            raise TypeError("'old_name' and 'new_name' must be strings.")
        if old_name not in points.headers:
            raise ValueError(f"Attribute '{old_name}' not found in dataset headers.")
        if old_name in ("id", "x", "y"):
            raise ValueError(f"Cannot rename core attribute '{old_name}'.")
        if new_name in ("id", "x", "y"):
            raise ValueError(f"Cannot rename to core attribute '{new_name}'.")
        if new_name in points.headers:
            raise ValueError(f"Attribute '{new_name}' already exists in dataset headers.")

        # Copy points (attrs unchanged, only headers matter here)
        new_points = [_Point(p.id, p.x, p.y, list(p.attrs)) for p in points]

        # Update headers
        new_headers = [new_name if h == old_name else h for h in points.headers]

        return _PointDataset(new_points, new_headers, id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    def attribute_as_list(self, points, attribute):
        """Return the specified attribute ('attribute') from each point in a PointDataset ('points') as a list. Any attribute including the point id or x or y coordinates can be returned.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered.
        
        attribute (string): the point attribute to return as a list."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if not isinstance(attribute, str):
            raise TypeError("'attribute' must be a string.")
        if attribute not in points.headers:
            raise ValueError(f"Attribute '{attribute}' not found in dataset headers.")

        values = []
        for p in points:
            if attribute == "id":
                values.append(p.id)
            elif attribute == "x":
                values.append(p.x)
            elif attribute == "y":
                values.append(p.y)
            else:
                # Offset by 3 since attrs start after id, x, y
                attr_index = points.headers.index(attribute) - 3
                values.append(p.attrs[attr_index])
        return values

    def spatial_join(self, points, to_join, columns="all",radius = None):
        """Join the data of one PointDataset ('points') with the data from another PointDataset ('to_join'), based on nearest neighbour spatial proximity. This join function is a left outer join, using a one-to-one nearest neighbour relationship only (i.e., other join types and relationships are not supported).
        
        Parameters:
        
        points (object): the PointDataset to join attributes onto. The input dataset will not be altered by the operation.
        
        to_join (object): the PointDataset to join from. Attributes from this dataset will be transferred.
        
        columns (list, string): a list of strings representing the attribute headers in 'to_join' to be joined to 'points'. May be specified as the string "all" to include all available attributes.
        
        radius (integer, float, None): optional parameter to only join points within a specified radius."""
        if not isinstance(points, _PointDataset) or not isinstance(to_join, _PointDataset):
            raise TypeError("Both 'points' and 'to_join' must be PointDataset objects.")

        # Normalize columns argument
        if isinstance(columns, str):
            columns = [columns]
        if columns[0] == "all":
            # Exclude id/x/y from joinable attributes
            columns = to_join.headers[:]

        # Validate requested columns
        for col_name in columns:
            if col_name not in to_join.headers:
                raise ValueError(f"Column '{col_name}' does not exist in to_join dataset.")

        # Copy points so we don't alter the original dataset
        new_points = [_Point(p.id, p.x, p.y, p.attrs.copy()) for p in points.points]

        # Extend headers with joined columns
        new_headers = points.headers.copy()
        for col_name in columns:
            if col_name not in new_headers:
                new_headers.append(col_name)

        for p in new_points:
            if radius is None:
                # Nearest neighbour join
                nearest = None
                nearest_dist = float("inf")
                for q in to_join.points:
                    dx = p.x - q.x
                    dy = p.y - q.y
                    dist = math.sqrt(dx * dx + dy * dy)
                    if dist < nearest_dist:
                        nearest = q
                        nearest_dist = dist
                if nearest is not None:
                    for col_name in columns:
                        attr_idx = to_join.headers.index(col_name) - 3
                        p.attrs.append(nearest.attrs[attr_idx])
                else:
                    for _ in columns:
                        p.attrs.append("")
            else:
                # Radius-based join
                joined_values = {col: [] for col in columns}
                for q in to_join.points:
                    dx = p.x - q.x
                    dy = p.y - q.y
                    dist = math.sqrt(dx * dx + dy * dy)
                    if dist <= radius:
                        for col_name in columns:
                            attr_idx = to_join.headers.index(col_name) - 3
                            joined_values[col_name].append(q.attrs[attr_idx])
                for col_name in columns:
                    p.attrs.append(joined_values[col_name] if joined_values[col_name] else "")

        return _PointDataset(new_points, new_headers, id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    # ============================================================
    # ATTRIBUTE FILTERING
    # ============================================================

    def filter_numeric_attribute(self, points, attribute, condition, threshold):
        """Remove points from a PointDataset ('points') that do not meet a condition placed upon a numerical attribute ('attribute') and return a new PointDataset.

        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        attribute (string): the point attribute to filter points by. Must be a header in the dataset for a numerical attribute.
        
        condition (string): the condition placed upon the point attribute for removing points. May be specified as "<" for less than, "<=" for less than equal to, ">" for greater than, ">=" for greater than equal to, "==" for equal to, or "!=" for not equal to.
        
        threshold (integer, float): the value to evaluate the condition against."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if condition not in ("<", "<=", ">", ">=", "==", "!="):
            raise ValueError("Condition must be one of '<', '<=', '>', '>=', '==', '!='.")
        
        ops = {
            "<": operator.lt,
            "<=": operator.le,
            ">": operator.gt,
            ">=": operator.ge,
            "==": operator.eq,
            "!=": operator.ne,
        }

        # Find the index of the attribute in headers
        try:
            attr_idx = points.headers.index(attribute)
        except ValueError:
            raise ValueError(f"Attribute '{attribute}' not found in dataset headers.")

        output_points = []
        for p in points.points:
            attr_val = p.attrs[attr_idx - 3]  # adjust index: headers include id/x/y first
            if ops[condition](attr_val, threshold):
                output_points.append(_Point(p.id, p.x, p.y, p.attrs.copy()))

        return _PointDataset(output_points, points.headers.copy(), id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    def filter_text_attribute(self, points, attribute, condition, threshold):
        """Filter a PointDataset ('points') by removing points that do not meet a specified threshold condition on a string (text) attribute, and return a new PointDataset.

        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        attribute (string): the point attribute to filter points by. Must be a header in the dataset for a string attribute.
        
        condition (string): the condition placed upon the point attribute for removing points. May be specified as "==" for equal to, "!=" for not equal to, "like" for contains, "not like" for does not contain, "start" for starts with, "not start" for does not start with, "end" for ends with, or "not end" for does not end with. Conditions are case sensitive.
        
        threshold (string): the value to evaluate the condition against."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")

        valid_conditions = {
            "==", "!=", "like", "not like", "start", "not start", "end", "not end"
        }
        if condition not in valid_conditions:
            raise ValueError(
                f"Invalid condition '{condition}'. Must be one of: {sorted(valid_conditions)}."
            )

        try:
            attr_idx = points.headers.index(attribute)
        except ValueError:
            raise ValueError(f"Attribute '{attribute}' not found in dataset headers.")

        output_points = []
        for p in points.points:
            val = p.attrs[attr_idx - 3]  # adjust index: headers include id/x/y first
            if not isinstance(val, str):
                continue  # non-string values fail the condition

            keep = False
            if condition == "==":
                keep = val == threshold
            elif condition == "!=":
                keep = val != threshold
            elif condition == "like":
                keep = threshold in val
            elif condition == "not like":
                keep = threshold not in val
            elif condition == "start":
                keep = val.startswith(threshold)
            elif condition == "not start":
                keep = not val.startswith(threshold)
            elif condition == "end":
                keep = val.endswith(threshold)
            elif condition == "not end":
                keep = not val.endswith(threshold)

            if keep:
                output_points.append(_Point(p.id, p.x, p.y, p.attrs.copy()))

        return _PointDataset(output_points, points.headers.copy(), id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    def filter_date_attribute(self, points, attribute, condition, threshold):
        """Filter a PointDataset ('points') by removing points that do not meet a specified threshold condition and return a new PointDataset. This filtering operation is performed with point attributes of date values, which must be in "yyyy-mm-dd" format.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        attribute (string): the point attribute to filter points by. Must be a header in the dataset for a date attribute in "yyyy-mm-dd" format.
        
        condition (string): the condition placed upon the point attribute for removing points. May be specified as "<" for less than, "<=" for less than equal to, ">" for greater than, ">=" for greater than equal to, "==" for equal to, or "!=" for not equal to.
        
        threshold (string): the threshold to evaluate the condition against. Must be in "yyyy-mm-dd" format."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")

        ops = {
            "<": operator.lt,
            "<=": operator.le,
            ">": operator.gt,
            ">=": operator.ge,
            "==": operator.eq,
            "!=": operator.ne,
        }

        if condition not in ops:
            raise ValueError(f"Invalid condition '{condition}'. Must be one of {list(ops.keys())}.")

        if not isinstance(threshold, str) or len(threshold) != 10:
            raise ValueError("Threshold must be a date string in 'YYYY-MM-DD' format.")

        try:
            attr_idx = points.headers.index(attribute)
        except ValueError:
            raise ValueError(f"Attribute '{attribute}' not found in dataset headers.")

        op = ops[condition]
        output_points = []
        for p in points.points:
            val = p.attrs[attr_idx - 3]  # adjust index: headers include id/x/y first
            if isinstance(val, str) and len(val) == 10:
                if op(val, threshold):
                    output_points.append(_Point(p.id, p.x, p.y, p.attrs.copy()))

        return _PointDataset(output_points, points.headers.copy(), id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    def filter_boolean_attribute(self, points, attribute, condition, target):
        """Remove points from a PointDataset ('points') that do not meet a Boolean condition placed upon an attribute ('attribute') and return a new PointDataset.

        Parameters:

        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.

        attribute (string): the point attribute to filter points by. Must be a header in the dataset for a Boolean attribute.

        condition (string): the condition placed upon the attribute for removing points. Must be "==" for equality or "!=" for inequality.

        target (bool): the Boolean value to evaluate the condition against."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")

        # Only equality and inequality make sense for booleans
        ops = {"==": operator.eq, "!=": operator.ne}

        if condition not in ops:
            raise ValueError(
                f"Invalid condition '{condition}'. Must be one of {list(ops.keys())}."
            )

        if not isinstance(target, bool):
            raise TypeError("Target must be a Boolean value (True or False).")

        # Find the index of the attribute in headers
        try:
            attr_idx = points.headers.index(attribute)
        except ValueError:
            raise ValueError(f"Attribute '{attribute}' not found in dataset headers.")

        output_points = []
        op = ops[condition]

        for p in points.points:
            attr_val = p.attrs[attr_idx - 3]  # adjust index: id/x/y come first

            # Only consider actual booleans
            if isinstance(attr_val, bool) and op(attr_val, target):
                output_points.append(_Point(p.id, p.x, p.y, p.attrs.copy()))

        return _PointDataset(
            output_points,
            points.headers.copy(),
            id_col=points.id_col,
            x_col=points.x_col,
            y_col=points.y_col
        )

    def filter_type_attribute(self, points, attribute, condition, dtype):
        """Filter a PointDataset ('points') by removing points that do not meet a specified data type condition and return a new PointDataset. This filtering operation is performed with point attributes of specific data types.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        attribute (string): the point attribute to filter points by. Must be a header in the dataset.
        
        condition (string): the condition placed upon the point attribute for removing points. May be specified as "is" to include values that match the data type, or "is not" to include values that do not match.
        
        dtype (string): the data type to evaluate each value in the specified attribute against to determine inclusion in the filtered data. Possible types include "int", "float", "string", "bool", and "nan"."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if condition not in ("is", "is not"):
            raise ValueError("Condition must be 'is' or 'is not'.")
        
        valid_types = {"int", "float", "string", "bool", "nan"}
        if dtype not in valid_types:
            raise ValueError(f"Invalid dtype '{dtype}'. Must be one of {valid_types}.")

        try:
            attr_idx = points.headers.index(attribute)
        except ValueError:
            raise ValueError(f"Attribute '{attribute}' not found in dataset headers.")

        def matches(val):
            if dtype == "int":
                return isinstance(val, int) and not isinstance(val, bool)  # exclude bools (since bool is subclass of int)
            elif dtype == "float":
                return isinstance(val, float) and not math.isnan(val)
            elif dtype == "string":
                return isinstance(val, str)
            elif dtype == "bool":
                return isinstance(val, bool)
            elif dtype == "nan":
                return isinstance(val, float) and math.isnan(val)
            return False

        output_points = []
        for p in points.points:
            val = p.attrs[attr_idx - 3]  # adjust index: headers include id/x/y first
            is_match = matches(val)
            keep = (is_match if condition == "is" else not is_match)
            if keep:
                output_points.append(_Point(p.id, p.x, p.y, p.attrs.copy()))

        return _PointDataset(output_points, points.headers.copy(), id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    def filter_missing_invalid_attributes(self, points, to_remove=""):
        """Filter a PointDataset ('points') by removing points that contain missing or invalid values ('to_remove') and return a new PointDataset.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered.
        
        to_remove (integer, float, string, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")

        # Normalize to list
        if not isinstance(to_remove, list):
            to_remove = [to_remove]

        output_points = []
        for p in points.points:
            # Skip points containing any invalid value in attrs
            if any(val in to_remove for val in p.attrs):
                continue
            output_points.append(_Point(p.id, p.x, p.y, p.attrs.copy()))

        return _PointDataset(output_points, points.headers.copy(), id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    # ============================================================
    # SPATIAL FILTERING AND EXTENTS
    # ============================================================

    def extent_area(self, points, decimals=2):
        """Calculate and return the area of the extent (bounding box) of a PointDataset ('points').
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        decimals (integer): the number of decimal places the area will be rounded to. Default is 2."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")

        min_x, min_y, max_x, max_y = points.extent
        area = (max_x - min_x) * (max_y - min_y)
        return round(abs(area), decimals)
    
    def extent_perimeter(self, points, decimals=2):
        """Calculate and return the perimeter of the extent (bounding box) of a PointDataset ('points').
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        decimals (integer): the number of decimal places the perimeter will be rounded to. Default is 2."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")

        min_x, min_y, max_x, max_y = points.extent
        perim = 2 * ((max_x - min_x) + (max_y - min_y))
        return round(abs(perim), decimals)

    def extent_center(self, points, decimals=2):
        """Calculate and return the center (midpoint) of the extent (bounding box) of a PointDataset ('points') as a tuple.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        decimals (integer): the number of decimal places the coordinates will be rounded to. Default is 2."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")

        min_x, min_y, max_x, max_y = points.extent
        cx = (min_x + max_x) / 2
        cy = (min_y + max_y) / 2
        return round(cx, decimals), round(cy, decimals)

    def point_distribution(self, points, iterations=100, decimals=2):
        """Calculate and return an index of the spatial distribution (pDist) of a PointDataset ('points').
        
        The pDist index compares the standard distance of the input dataset to the average standard distance
        of random point distributions within the same extent. Values < 0 indicate clustering, ~0 randomness,
        and > 0 dispersion.
        
        Parameters:
        
        points (object): the dataset to analyze. The input dataset will not be altered by the operation.
        
        iterations (integer): number of random point distributions to generate for reference. Default is 100.
        
        decimals (integer): number of decimal places to round values. Default is 2."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if len(points) == 0:
            raise ValueError("Cannot calculate distribution of an empty dataset.")
        if not isinstance(iterations, int) or iterations <= 0:
            raise ValueError("'iterations' must be a positive integer.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")

        # observed standard distance
        target_mean = self.standard_distance(points, decimals)
        np = len(points)

        # extent from dataset
        min_x, min_y, max_x, max_y = points.extent

        # Monte Carlo reference
        test_dist_sum = 0
        for _ in range(iterations):
            test_points = TableTools().generate.rand_points((min_x, min_y, max_x, max_y), np, decimals)
            test_dist = self.standard_distance(test_points, decimals)
            test_dist_sum += test_dist

        avg_test_mean = round(test_dist_sum / iterations, decimals)

        # pDist index
        pDist = round((target_mean - avg_test_mean) / (target_mean + avg_test_mean), decimals)
        return pDist
        
    def point_in_extent(self, point, extent):
        """Determine if a specified Point object ('point') falls within the specified bounding box extent ('extent').
        
        Parameters:
        
        point (object): the Point object to test.
        
        extent (list, tuple): a list or tuple representing the extent in the format (min_x, min_y, max_x, max_y)."""
        if not isinstance(extent, (list,tuple)) or len(extent) != 4:
            raise ValueError("'extent' must be a list or tuple of length 4 (min_x, min_y, max_x, max_y).")
        
        if not isinstance(point,_Point):
            raise ValueError("'point' must be a Point object")
        
        min_x, min_y, max_x, max_y = extent
        return (min_x <= point.x <= max_x) and (min_y <= point.y <= max_y)

    def clip_points_by_extent(self, points, extent, method):
        """Remove points from a PointDataset ('points') based on a specified output bounding box extent ('extent') and return a new PointDataset.

        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        extent (list): a list representing the minimum x, minimum y, maximum x, and maximum y, in that order, of the extent to remove points from.
        
        method (string): the method of removing points from the extent. May be specified as "inside" to remove points that fall inside of the specified extent, or "outside" to remove points that fall outside of the specified extent."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if not (isinstance(extent, (list, tuple)) and len(extent) == 4):
            raise ValueError("'extent' must be a list or tuple of four values [min_x, min_y, max_x, max_y].")
        if method not in ("inside", "outside"):
            raise ValueError("Method must be 'inside' or 'outside'.")

        minx, miny, maxx, maxy = extent

        output_points = []
        for p in points.points:
            px, py = p.x, p.y
            if method == "outside":
                # keep points inside the extent
                if minx < px < maxx and miny < py < maxy:
                    output_points.append(_Point(p.id, p.x, p.y, p.attrs.copy()))
            elif method == "inside":
                # keep points outside the extent
                if px < minx or px > maxx or py < miny or py > maxy:
                    output_points.append(_Point(p.id, p.x, p.y, p.attrs.copy()))

        return _PointDataset(output_points, points.headers.copy(), id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    # ============================================================
    # COORDINATE NORMALIZATION AND CONVERSION
    # ============================================================

    def convert_dms_to_dd(self, points, d_sep="°", m_sep="'", s_sep='"'):
        """Convert each coordinate value in a PointDataset ('points') from degrees, minutes, seconds to decimal degrees and return a new PointDataset.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        d_sep (string): the character that separates the degree component of the coordinate from the other components. The default value is °.
        
        m_sep (string): the character that separates the minute component of the coordinate from the other components. The default value is '.
        
        s_sep (string): the character that separates the second component of the coordinate from the other components. The default value is "."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")

        def dms_to_dd(coord: str) -> float:
            if not isinstance(coord, str):
                raise TypeError("Coordinate must be a string in DMS format.")
            # Normalize separators
            c = coord.replace(d_sep, ".").replace(m_sep, ".").replace(s_sep, "").replace(" ", "")
            parts = c.split(".")
            if len(parts) < 3:
                raise ValueError(f"Invalid DMS coordinate format: '{coord}'")
            deg, minutes, seconds = parts[0], parts[1], parts[2]
            return round(float(deg) + float(minutes)/60 + float(seconds)/3600, 6)

        # Convert each point
        new_points = []
        for p in points:
            new_x = dms_to_dd(p.x)
            new_y = dms_to_dd(p.y)
            new_points.append(_Point(p.id, new_x, new_y, list(p.attrs)))

        # Headers remain unchanged
        new_headers = list(points.headers)

        return _PointDataset(new_points, new_headers, id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)
    
    def convert_dd_to_dms(self, points, d_sep="°", m_sep="'", s_sep='"'):
        """Convert each coordinate value in a PointDataset ('points') from decimal degrees to degrees, minutes, seconds and return a new PointDataset.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        d_sep (string): the character that will separate the degree component of the coordinate from the other components. The default value is °.
        
        m_sep (string): the character that will separate the minute component of the coordinate from the other components. The default value is '.
        
        s_sep (string): the character that will separate the second component of the coordinate from the other components. The default value is "."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")

        def dd_to_dms(coord: float) -> str:
            if not isinstance(coord, (int, float)):
                raise TypeError("Coordinate must be numeric (decimal degrees).")
            d = int(coord)
            m = int((coord - d) * 60)
            s = round((coord - d - m/60) * 3600)

            d_str = f"{d:02d}{d_sep}"
            m_str = f"{m:02d}{m_sep}"
            s_str = f"{s:02d}{s_sep}"
            return d_str + m_str + s_str

        # Convert each point
        new_points = []
        for p in points:
            new_x = dd_to_dms(p.x)
            new_y = dd_to_dms(p.y)
            new_points.append(_Point(p.id, new_x, new_y, list(p.attrs)))

        # Headers remain unchanged
        new_headers = list(points.headers)

        return _PointDataset(new_points, new_headers, id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    def normalize_latitude(self, points):
        """Normalize latitude (y) values in a PointDataset ('points') to the range [-90, 90].
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")

        new_points = []
        for p in points:
            lat = p.y
            # Clamp to valid latitude range
            if lat > 90:
                lat = 90
            elif lat < -90:
                lat = -90
            new_points.append(_Point(p.id, p.x, lat, list(p.attrs)))

        return _PointDataset(new_points, list(points.headers), id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    def normalize_longitude(self, points, mode="180"):
        """Normalize longitude (x) values in a PointDataset ('points') to either [-180, 180] or [0, 360].
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered.
        
        mode (string): '180' for range [-180, 180], '360' for range [0, 360]. Default is '180'."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if mode not in ("180", "360"):
            raise ValueError("Mode must be '180' or '360'.")

        new_points = []
        for p in points:
            lon = p.x
            if mode == "180":
                lon = ((lon + 180) % 360) - 180
            else:  # mode == "360"
                lon = lon % 360
            new_points.append(_Point(p.id, lon, p.y, list(p.attrs)))

        return _PointDataset(new_points, list(points.headers), id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    # ============================================================
    # COORDINATE TRANSFORMATIONS
    # ============================================================

    def unrotate_x(self, x, y, pole_x, pole_y):
        """Unrotate a list of x coordinates ('x') from a rotated pole coordinate system and return a new list of x coordinates.
        
        Parameters:
        
        x (list): the list of x coordinates the operation will be performed on. The input list will not be altered by the operation.
        
        y (list): the list of y coordinates required for the operation.
        
        pole_x (integer, float): the x coordinate of the pole about which the x will be rotated.
        
        pole_y (integer, float): the y coordinate of the pole about which the x will be rotated."""
        if not isinstance(x, list) or not isinstance(y, list):
            raise TypeError("'x' and 'y' must be lists of numeric values.")
        if len(x) != len(y):
            raise ValueError("'x' and 'y' must be the same length.")

        rx = [math.radians(v) for v in x]
        ry = [math.radians(v) for v in y]

        s1 = math.sin(math.radians(pole_y))
        c1 = math.cos(math.radians(pole_y))
        s2 = math.sin(math.radians(pole_x))
        c2 = math.cos(math.radians(pole_x))

        unrotated_x = []
        for xv, yv in zip(rx, ry):
            tmp1 = s2 * (-s1 * math.cos(xv) * math.cos(yv) + c1 * math.sin(yv)) - c2 * math.sin(xv) * math.cos(yv)
            tmp2 = c2 * (-s1 * math.cos(xv) * math.cos(yv) + c1 * math.sin(yv)) + s2 * math.sin(xv) * math.cos(yv)

            # atan2 handles quadrant correctly and avoids division by zero
            final_x = math.degrees(math.atan2(tmp1, tmp2))
            unrotated_x.append(final_x)

        return unrotated_x

    def unrotate_y(self, x, y, pole_x, pole_y):
        """Unrotate a list of y coordinates ('y') from a rotated pole coordinate system and return a new list of y coordinates.
        
        Parameters:
        
        x (list): the list of x coordinates required for the operation.
        
        y (list): the list of y coordinates the operation will be performed on. The input list will not be altered by the operation.
        
        pole_x (integer, float): the x coordinate of the pole about which the y will be rotated.
        
        pole_y (integer, float): the y coordinate of the pole about which the y will be rotated."""
        if not isinstance(x, list) or not isinstance(y, list):
            raise TypeError("'x' and 'y' must be lists of numeric values.")
        if len(x) != len(y):
            raise ValueError("'x' and 'y' must be the same length.")

        rx = [math.radians(v) for v in x]
        ry = [math.radians(v) for v in y]

        s1 = math.sin(math.radians(pole_y))
        c1 = math.cos(math.radians(pole_y))

        unrotated_y = []
        for xv, yv in zip(rx, ry):
            val = s1 * math.sin(yv) + c1 * math.cos(yv) * math.cos(xv)

            # Clamp to [-1, 1] before asin to avoid domain errors due to floating-point drift
            val = max(-1.0, min(1.0, val))

            final_y = math.degrees(math.asin(val))
            unrotated_y.append(final_y)

        return unrotated_y

    def unrotate_coordinates(self, points, pole_x, pole_y):
        """Unrotate the x and y coordinates in a PointDataset ('points') from a rotated pole coordinate system and return a new PointDataset. The x and y coordinates will be overwritten with the unrotated coordinates.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        pole_x (integer, float): the x coordinate of the pole about which the points will be rotated.
        
        pole_y (integer, float): the y coordinate of the pole about which the points will be rotated."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")

        # Extract coordinate lists
        xs = [p.x for p in points]
        ys = [p.y for p in points]

        # Apply unrotation
        xr = self.unrotate_x(xs, ys, pole_x, pole_y)
        yr = self.unrotate_y(xs, ys, pole_x, pole_y)

        # Build new points with unrotated coordinates
        new_points = []
        for i, p in enumerate(points):
            new_points.append(_Point(p.id, xr[i], yr[i], list(p.attrs)))

        # Headers remain unchanged
        new_headers = list(points.headers)

        return _PointDataset(new_points, new_headers, id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    def translate(self, points, x_shift=0.0, y_shift=0.0):
        """Displace the x and y coordinates of a PointDataset ('points') by a specified x and y offset distance and return a new PointDataset.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        x_shift (integer, float): the distance to displace the x coordinate of each point.
        
        y_shift (integer, float): the distance to displace the y coordinate of each point."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if not isinstance(x_shift, (int, float)) or not isinstance(y_shift, (int, float)):
            raise TypeError("'x_shift' and 'y_shift' must be numeric.")

        new_points = []
        for p in points:
            new_points.append(_Point(p.id, p.x + x_shift, p.y + y_shift, list(p.attrs)))

        return _PointDataset(new_points, list(points.headers), id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    def translate_random_uniform(self, points, dist):
        """Randomly displace the x and y coordinates of a PointDataset ('points') by a uniform random offset between -dist and +dist, and return a new PointDataset.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        dist (integer, float): the maximum absolute distance to randomly displace the coordinates of each point."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if not isinstance(dist, (int, float)):
            raise TypeError("'dist' must be numeric.")

        new_points = []
        for p in points:
            dx = random.uniform(-dist, dist)
            dy = random.uniform(-dist, dist)
            new_points.append(_Point(p.id, p.x + dx, p.y + dy, list(p.attrs)))

        return _PointDataset(new_points, list(points.headers), id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    def scale(self, points, x_factor=1.0, y_factor=1.0):
        """Scale the x and y coordinates of a PointDataset ('points') by specified factors and return a new PointDataset.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation. 
        
        x_factor (integer, float): the factor to multiply the x coordinate of each point. Default is 1.0.
        
        y_factor (integer, float): the factor to multiply the y coordinate of each point. Default is 1.0."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if not isinstance(x_factor, (int, float)) or not isinstance(y_factor, (int, float)):
            raise TypeError("'x_factor' and 'y_factor' must be numeric.")

        new_points = []
        for p in points:
            new_points.append(_Point(p.id, p.x * x_factor, p.y * y_factor, list(p.attrs)))

        return _PointDataset(new_points, list(points.headers), id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    def rotate(self, points, angle, pivot_x=0.0, pivot_y=0.0):
        """Rotate all points in a PointDataset ('points') around a specified pivot by a given angle (in degrees) and return a new PointDataset.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        angle (integer, float): the rotation angle in degrees. Positive values rotate counterclockwise.
        
        pivot_x (integer, float): the x coordinate of the pivot point. Default is 0.0.
        
        pivot_y (integer, float): the y coordinate of the pivot point. Default is 0.0."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if not isinstance(angle, (int, float)):
            raise TypeError("'angle' must be numeric.")
        if not isinstance(pivot_x, (int, float)) or not isinstance(pivot_y, (int, float)):
            raise TypeError("'pivot_x' and 'pivot_y' must be numeric.")

        # Convert angle to radians
        theta = math.radians(angle)
        cos_t = math.cos(theta)
        sin_t = math.sin(theta)

        new_points = []
        for p in points:
            # Translate point relative to pivot
            dx = p.x - pivot_x
            dy = p.y - pivot_y

            # Apply rotation
            new_x = dx * cos_t - dy * sin_t + pivot_x
            new_y = dx * sin_t + dy * cos_t + pivot_y

            new_points.append(_Point(p.id, new_x, new_y, list(p.attrs)))

        return _PointDataset(new_points, list(points.headers), id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    def reflect(self, points, axis="x"):
        """Reflect all points in a PointDataset ('points') across a specified axis and return a new PointDataset.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation..
        
        axis (string): the axis across which to reflect the points. Options are 'x' (mirror across x-axis), 'y' (mirror across y-axis), or 'origin' (mirror across both axes). Default is 'x'."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if axis not in ("x", "y", "origin"):
            raise ValueError("Axis must be 'x', 'y', or 'origin'.")

        new_points = []
        for p in points:
            if axis == "x":
                new_points.append(_Point(p.id, p.x, -p.y, list(p.attrs)))
            elif axis == "y":
                new_points.append(_Point(p.id, -p.x, p.y, list(p.attrs)))
            elif axis == "origin":
                new_points.append(_Point(p.id, -p.x, -p.y, list(p.attrs)))

        return _PointDataset(new_points, list(points.headers), id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    def shear(self, points, x_shear=0.0, y_shear=0.0):
        """Apply a shear transformation to all points in a PointDataset ('points') and return a new PointDataset.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        x_shear (integer, float): shear factor applied in the x direction (skews x by a proportion of y). Default is 0.0.
        
        y_shear (integer, float): shear factor applied in the y direction (skews y by a proportion of x). Default is 0.0."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if not isinstance(x_shear, (int, float)) or not isinstance(y_shear, (int, float)):
            raise TypeError("'x_shear' and 'y_shear' must be numeric.")

        new_points = []
        for p in points:
            new_x = p.x + x_shear * p.y
            new_y = p.y + y_shear * p.x
            new_points.append(_Point(p.id, new_x, new_y, list(p.attrs)))

        return _PointDataset(new_points, list(points.headers), id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    def affine_transformation(self, points, dx=0.0, dy=0.0, x_factor=1.0, y_factor=1.0, angle=0.0, pivot_x=0.0, pivot_y=0.0, reflect_axis=None, x_shear=0.0, y_shear=0.0):
        """Apply a general affine transformation (translate, scale, rotate, reflect, shear) to a PointDataset ('points') and return a new PointDataset.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on.  The input dataset will not be altered by the operation.
        
        dx (float): translation offsets in the x direction. Default is 0.0.

        dy (float): translation offset in the y direction. Default is 0.0.
        
        x_factor (float): scaling factor for x coordinate. Default is 1.0.

        y_factor (float): scaling factor for y coordinate. Default is 1.0.
        
        angle (float): rotation angle in degrees (counterclockwise). Default is 0.0.
        
        pivot_x (float): pivot point for rotation. Default is 0.0,.

        pivot_y (float): pivot point for rotation. Default is 0.0.
        
        reflect_axis (string): optional reflection axis ('x', 'y', 'origin'). Default is None.
        
        x_shear (float): shear factor in x direction. Default is 0.0.

        y_shear (float): shear factor in x y directions. Default is 0.0."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")

        theta = math.radians(angle)
        cos_t, sin_t = math.cos(theta), math.sin(theta)

        new_points = []
        for p in points:
            x, y = p.x, p.y

            # 1. Translate
            x += dx
            y += dy

            # 2. Scale
            x *= x_factor
            y *= y_factor

            # 3. Rotate (around pivot)
            dx_rel, dy_rel = x - pivot_x, y - pivot_y
            x = dx_rel * cos_t - dy_rel * sin_t + pivot_x
            y = dx_rel * sin_t + dy_rel * cos_t + pivot_y

            # 4. Reflect
            if reflect_axis == "x":
                y = -y
            elif reflect_axis == "y":
                x = -x
            elif reflect_axis == "origin":
                x, y = -x, -y

            # 5. Shear
            x = x + x_shear * y
            y = y + y_shear * x

            new_points.append(_Point(p.id, x, y, list(p.attrs)))

        return _PointDataset(new_points, list(points.headers), id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)
   
    def swap_coordinates(self, points):
        """Switch the x and y coordinates in a PointDataset ('points') and return a new PointDataset.
        
        Parameters:
        
        points (points): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")

        # Copy points with swapped coordinates
        new_points = []
        for p in points:
            new_points.append(_Point(p.id, p.y, p.x, list(p.attrs)))

        # Headers remain the same (id, x, y, plus attrs)
        new_headers = list(points.headers)

        return _PointDataset(new_points, new_headers, id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    # ============================================================
    # DISTANCE AND NEIGHBOURHOOD OPERATIONS
    # ============================================================

    def nearest_neighbour(self, points1, points2, num_nearest=1, method="Euclidean", decimals=2):
        """Compute the nearest neighbours of each point in an input PointDataset ('points1') from a second PointDataset ('points2') and return a new PointDataset with added headers and attributes for nearest neighbour IDs and distances.

        Parameters:
        
        points1 (object): the dataset to determine nearest neighbours of. The input dataset will not be altered by the operation.
        
        points2 (object): the dataset to determine nearest neighbours from.
        
        num_nearest (integer): the number of nearest neighbours to determine for each point.
        
        method (string): the method with which to calculate distance. Possible choices are 'Euclidean' and 'Manhattan'.
        
        decimals (integer): the number of decimal places distance will be rounded to. Default is 2."""
        if not isinstance(points1, _PointDataset) or not isinstance(points2, _PointDataset):
            raise TypeError("'points1' and 'points2' must be PointDataset objects.")
        if not isinstance(num_nearest, int) or num_nearest <= 0:
            raise ValueError("'num_nearest' must be a positive integer.")
        if method not in ("Euclidean", "Manhattan"):
            raise ValueError("Method must be 'Euclidean' or 'Manhattan'.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")

        # Copy points1 so we don't alter the original dataset
        new_points = [_Point(p.id, p.x, p.y, p.attrs.copy()) for p in points1.points]

        # First pass: compute neighbours per point without changing headers
        neighbours_per_point = []
        global_k = 0  # maximum available neighbours across all points

        same_dataset = (points1 is points2)

        for p in new_points:
            results = []
            for p2 in points2.points:
                # avoid self-match only when datasets are identical
                if not (same_dataset and p.id == p2.id):
                    if method == "Euclidean":
                        dx = p.x - p2.x
                        dy = p.y - p2.y
                        dist = round(math.sqrt(dx * dx + dy * dy), decimals)
                    else:  # Manhattan
                        dx = abs(p.x - p2.x)
                        dy = abs(p.y - p2.y)
                        dist = round(dx + dy, decimals)
                    results.append((p2.id, dist))

            results.sort(key=lambda x: x[1])
            k_i = min(num_nearest, len(results))
            global_k = max(global_k, k_i)
            neighbours_per_point.append(results[:k_i])

        # If no neighbours exist for any point, return dataset unchanged
        if global_k == 0:
            return _PointDataset(new_points, points1.headers.copy(),id_col=points1.id_col, x_col=points1.x_col, y_col=points1.y_col)

        # Second pass: add headers and populate attrs
        new_headers = points1.headers.copy()
        for n in range(1, global_k + 1):
            new_headers.append(f'nearest_id_{n}')
            new_headers.append(f'nearest_dist_{n}')

        # Append values in the same order as headers; pad with None if a point has fewer neighbours
        for p, neigh in zip(new_points, neighbours_per_point):
            # add available neighbours
            for n in range(len(neigh)):
                p.attrs.append(neigh[n][0])  # nearest_id_{n+1}
                p.attrs.append(neigh[n][1])  # nearest_dist_{n+1}
            # pad missing neighbours up to global_k
            for _ in range(len(neigh), global_k):
                p.attrs.append(None)  # nearest_id_{missing}
                p.attrs.append(None)  # nearest_dist_{missing}

        return _PointDataset(new_points, new_headers, id_col=points1.id_col, x_col=points1.x_col, y_col=points1.y_col)
 
    def distance(self, points1, points2, method="Euclidean", decimals=2):
        """Compute the distance between each point in an input PointDataset ('points1') and the corresponding point in a second PointDataset ('points2'), and return a new PointDataset with an added 'dist' column.

        Parameters:
        
        points1 (object): the dataset to calculate distances for. The input dataset will not be altered by the operation.
        
        points2 (object): the dataset of corresponding points to calculate distances to.
        
        method (string): the method with which to calculate distance. Possible choices are 'Euclidean' and 'Manhattan'.
        
        decimals (integer): the number of decimal places float values will be rounded to. Default is 2."""
        if not isinstance(points1, _PointDataset) or not isinstance(points2, _PointDataset):
            raise TypeError("'points1' and 'points2' must be PointDataset objects.")
        if len(points1) != len(points2):
            raise ValueError("Both datasets must have the same number of points.")
        if method not in ("Euclidean", "Manhattan"):
            raise ValueError("Method must be 'Euclidean' or 'Manhattan'.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")

        # Copy points1 so we don't alter the original dataset
        new_points = [_Point(p.id, p.x, p.y, p.attrs.copy()) for p in points1.points]

        # Extend headers with 'dist'
        new_headers = points1.headers.copy()
        new_headers.append("dist")

        # Compute distances
        for i, p in enumerate(new_points):
            p2 = points2.points[i]
            if method == "Euclidean":
                dx = p.x - p2.x
                dy = p.y - p2.y
                dist = round(math.sqrt(dx * dx + dy * dy), decimals)
            else:  # Manhattan
                dx = abs(p.x - p2.x)
                dy = abs(p.y - p2.y)
                dist = round(dx + dy, decimals)
            p.attrs.append(dist)

        return _PointDataset(new_points, new_headers, id_col=points1.id_col, x_col=points1.x_col, y_col=points1.y_col)

    def num_points_in_range(self, points, radius):
        """Calculate the number of points around each point in a PointDataset ('points') within a specified distance ('radius') and return a new PointDataset.

        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        radius (integer, float): the distance to search around each point for other points."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if not isinstance(radius, (int, float)):
            raise ValueError("'radius' must be a number.")

        # Copy points so we don't alter the original dataset
        new_points = [_Point(p.id, p.x, p.y, p.attrs.copy()) for p in points.points]

        # Extend headers with 'in_range'
        new_headers = points.headers.copy()
        new_headers.append("in_range")

        # Compute number of points in range for each point
        for i, p1 in enumerate(new_points):
            num_points = 0
            for j, p2 in enumerate(new_points):
                if i != j:  # avoid self-comparison
                    dx = p1.x - p2.x
                    dy = p1.y - p2.y
                    dist = math.sqrt(dx * dx + dy * dy)
                    if dist <= radius:
                        num_points += 1
            p1.attrs.append(num_points)

        return _PointDataset(new_points, new_headers, id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    # ============================================================
    # INTERPOLATION
    # ============================================================

    def interpolate_inverse_distance_weighted(self, points, attr_name, radius, power=2, to_fill="", decimals=2):
        """Interpolate a new value for each missing or invalid value ('to_fill') for a specified attribute ('attr_name') in a PointDataset ('points') using an inverse distance weighted (IDW) technique and return a new PointDataset.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        attr_name (string): the name of the attribute to interpolate. Only missing or invalid values will be interpolated.
        
        radius (integer, float): the distance to search around each point for other points to interpolate with. Only valid (numeric) values of points within this range will be used for interpolation.
        
        power (integer): the power of the inverse distance weighting function. Higher values will give greater weight to nearer values i.e. nearer values will contribute more strongly to the interpolated value.
        
        to_fill (integer, float, string, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value.
        
        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if not isinstance(radius, (int, float)):
            raise ValueError("'radius' must be a number.")
        if not isinstance(power, int) or power <= 0:
            raise ValueError("'power' must be a positive integer.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")

        if not isinstance(to_fill, list):
            to_fill = [to_fill]

        try:
            attr_idx = points.headers.index(attr_name)
        except ValueError:
            raise ValueError(f"Attribute '{attr_name}' not found in dataset headers.")

        # Copy points so we don't alter the original dataset
        new_points = [_Point(p.id, p.x, p.y, p.attrs.copy()) for p in points.points]

        # IDW interpolation
        for i, p1 in enumerate(new_points):
            val = p1.attrs[attr_idx - 3]  # adjust index: headers include id/x/y first
            if val in to_fill:  # if we need to interpolate this value
                numerators = []
                denominators = []
                for j, p2 in enumerate(new_points):
                    if i != j:  # don't use the same point
                        val2 = p2.attrs[attr_idx - 3]
                        if val2 not in to_fill and isinstance(val2, (int, float)):
                            dx = p1.x - p2.x
                            dy = p1.y - p2.y
                            dist = math.sqrt(dx * dx + dy * dy)
                            if dist <= radius and dist > 0:
                                numerators.append(val2 / (dist ** power))
                                denominators.append(1 / (dist ** power))
                if numerators and denominators:
                    p1.attrs[attr_idx - 3] = round(sum(numerators) / sum(denominators), decimals)

        return _PointDataset(new_points, points.headers.copy(),
                            id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    # ============================================================
    # CLUSTERING
    # ============================================================

    def kmeans_cluster(self, points, k, attributes=["x", "y"], convergence_threshold=1e-5, max_iterations=50):
        """Cluster points in a PointDataset ('points') into k clusters using the k-means algorithm and return a new PointDataset with a 'cluster' attribute appended.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        k (integer): the number of clusters to generate. Must be less than or equal to the number of points.
        
        attributes (list): list of attribute headers to use for clustering. Must be numeric attributes. Must include at least two attributes.
        
        convergence_threshold (float): the threshold for centroid movement to determine convergence.
        
        max_iterations (integer): the maximum number of iterations to perform. """
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if not isinstance(k, int) or k <= 1:
            raise ValueError("'k' must be an integer greater than 1.")
        if k > len(points.points):
            raise ValueError("'k' cannot be greater than the number of points in the dataset.")
        if not isinstance(attributes, list) or len(attributes) < 2:
            raise ValueError("'attributes' must be a list of at least two headers.")

        # Validate attributes exist
        for attr in attributes:
            if attr not in points.headers:
                raise ValueError(f"Attribute '{attr}' not found in dataset headers.")

        # Helper: extract attribute vector from a point
        def get_vector(p):
            vec = []
            for attr in attributes:
                if attr == points.x_col:
                    vec.append(p.x)
                elif attr == points.y_col:
                    vec.append(p.y)
                elif attr == points.id_col:
                    raise ValueError("ID column cannot be used for clustering.")
                else:
                    idx = points.headers.index(attr) - 3
                    val = p.attrs[idx]
                    if not isinstance(val, (int, float)):
                        raise ValueError(f"Attribute '{attr}' must be numeric.")
                    vec.append(val)
            return vec

        # Initialize centroids randomly
        centroids = [get_vector(p) for p in random.sample(points.points, k)]

        # Copy points so we don't alter the original dataset
        new_points = [_Point(p.id, p.x, p.y, p.attrs.copy()) for p in points.points]

        # Extend headers with 'cluster'
        new_headers = points.headers.copy()
        if "cluster" not in new_headers:
            new_headers.append("cluster")

        converged = False
        iteration = 0

        while not converged and iteration < max_iterations:
            iteration += 1

            # Assign points to nearest centroid
            for p in new_points:
                vec = get_vector(p)
                distances = [math.sqrt(sum((vx - cx) ** 2 for vx, cx in zip(vec, centroid))) for centroid in centroids]
                cluster_id = distances.index(min(distances))
                # Ensure cluster label is last attr
                if len(p.attrs) < len(new_headers) - 3:
                    p.attrs.append(cluster_id)
                else:
                    p.attrs[-1] = cluster_id

            # Recompute centroids
            new_centroids = []
            total_shift = 0.0
            for cluster_id in range(k):
                cluster_points = [get_vector(p) for p in new_points
                                if p.attrs[-1] == cluster_id]
                if cluster_points:
                    # Mean of cluster points
                    dim_means = [sum(dim) / len(cluster_points)
                                for dim in zip(*cluster_points)]
                    new_centroids.append(dim_means)
                    # Compute shift
                    dx = sum((a - b) ** 2 for a, b in zip(dim_means, centroids[cluster_id]))
                    total_shift += math.sqrt(dx)
                else:
                    # Empty cluster, reinitialize randomly
                    new_centroids.append(get_vector(random.choice(new_points)))

            centroids = new_centroids
            converged = total_shift < convergence_threshold

        return _PointDataset(new_points, new_headers, id_col=points.id_col, x_col=points.x_col, y_col=points.y_col)

    def dbscan_cluster(self, points, eps, min_pts, attributes=["x", "y"]):
        """Cluster points in a PointDataset ('points') using the DBSCAN algorithm and return a new PointDataset with a 'cluster' attribute appended. Noise points are assigned cluster ID -1.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        eps (float): the maximum radius distance to consider points as neighbours.
        
        min_pts (integer): the minimum number of points required within 'eps' distance for a point to be considered a core point.
        
        attributes (list): list of attribute headers to use for clustering. Must be numeric attributes. Must include at least two attributes."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if not isinstance(eps, (int, float)) or eps <= 0:
            raise ValueError("'eps' must be a positive number.")
        if not isinstance(min_pts, int) or min_pts < 2:
            raise ValueError("'min_pts' must be an integer greater than or equal to 2.")
        if min_pts > len(points.points):
            raise ValueError("'min_pts' cannot be greater than the number of points in the dataset.")
        if not isinstance(attributes, list) or len(attributes) < 2:
            raise ValueError("'attributes' must be a list of at least two headers.")

        # Validate attributes exist
        for attr in attributes:
            if attr not in points.headers:
                raise ValueError(f"Attribute '{attr}' not found in dataset headers.")

        # Helper: extract attribute vector from a point
        def get_vector(p):
            vec = []
            for attr in attributes:
                if attr == points.x_col:
                    vec.append(p.x)
                elif attr == points.y_col:
                    vec.append(p.y)
                elif attr == points.id_col:
                    raise ValueError("ID column cannot be used for clustering.")
                else:
                    idx = points.headers.index(attr) - 3
                    val = p.attrs[idx]
                    if not isinstance(val, (int, float)):
                        raise ValueError(f"Attribute '{attr}' must be numeric.")
                    vec.append(val)
            return vec

        # Helper: compute Euclidean distance between two vectors
        def dist(v1, v2):
            return math.sqrt(sum((a - b) ** 2 for a, b in zip(v1, v2)))

        # Helper: find neighbours of a point within eps
        def get_neighbours(target, all_points):
            target_vec = get_vector(target)
            neighbours = []
            for q in all_points:
                if q is not target:
                    if dist(target_vec, get_vector(q)) <= eps:
                        neighbours.append(q)
            return neighbours
        
        # Copy points so we don't alter the original dataset
        new_points = [_Point(p.id, p.x, p.y, p.attrs.copy()) for p in points.points]

        # Extend headers with 'cluster'
        new_headers = points.headers.copy()
        if "cluster" not in new_headers:
            new_headers.append("cluster")

        cluster_id = 0

        # Track visited points
        visited = set()

        for p in new_points:
            if p in visited:
                continue
            visited.add(p)

            neighbours = get_neighbours(p, new_points)

            if len(neighbours) < min_pts:
                # Not enough neighbours = noise
                p.attrs.append(-1)
            else:
                # Start a new cluster
                p.attrs.append(cluster_id)

                # Expand cluster
                queue = neighbours[:]
                while queue:
                    q = queue.pop(0)
                    if q not in visited:
                        visited.add(q)
                        q_neighbours = get_neighbours(q, new_points)
                        if len(q_neighbours) >= min_pts:
                            queue.extend(q_neighbours)

                    # Assign cluster if not already assigned
                    if len(q.attrs) < len(new_headers) - 3:
                        q.attrs.append(cluster_id)
                    elif q.attrs[-1] == -1:
                        q.attrs[-1] = cluster_id

                cluster_id += 1
        # Ensure every point has a cluster label
        for p in new_points:
            if len(p.attrs) < len(new_headers) - 3:
                p.attrs.append(-1)  # default to noise if never assigned

        # Return new dataset with cluster assignments
        return _PointDataset(new_points, new_headers,id_col=points.id_col,x_col=points.x_col,y_col=points.y_col)

    def cluster_summary(self, points):
        """Return a dictionary of cluster IDs and their counts from a clustered PointDataset.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. Must contain a 'cluster' attribute."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if "cluster" not in points.headers:
            raise ValueError("Dataset does not contain a 'cluster' attribute.")
        
        cluster_idx = points.headers.index("cluster") - 3  # adjust for id/x/y offset

        counts = {}
        for p in points.points:
            cid = p.attrs[cluster_idx]  # cluster label is last attribute
            counts[cid] = counts.get(cid, 0) + 1

        return dict(sorted(counts.items(), key=lambda x: x[0]))

    def medoid(self, points):
        """Compute and return the medoid (a central point that must belong to the set of input points) of a PointDataset ('points') as a Point object.
        
        Parameters:
        
        points (PointDataset): the dataset to analyze. The input dataset will not be altered by the operation."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if len(points) == 0:
            raise ValueError("Cannot compute medoid of an empty dataset.")

        # calculate the median coordinates
        xs = [p.x for p in points]
        ys = [p.y for p in points]
        xmed = statistics.median(xs)
        ymed = statistics.median(ys)

        # find the nearest neighbour of the median point
        nearest = []
        for p in points:
            dx = xmed - p.x
            dy = ymed - p.y
            dist = math.sqrt(dx * dx + dy * dy)
            nearest.append((p, dist))

        # sort by distance and return the closest point
        nearest.sort(key=lambda x: x[1])
        return nearest[0][0]

    # ============================================================
    # SPATIAL STATISTICS
    # ============================================================
          
    def mean_center(self, points, decimals=2):
        """Compute and return the mean center (average x and y coordinates) of a PointDataset ('points') as a tuple.
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        decimals (integer): the number of decimal places output coordinates will be rounded to. Default is 2."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")
        if len(points) == 0:
            raise ValueError("Cannot compute mean center of an empty dataset.")

        x_sum = sum(p.x for p in points)
        y_sum = sum(p.y for p in points)
        lp = len(points)

        return round(x_sum / lp, decimals), round(y_sum / lp, decimals)
      
    def weighted_mean_center(self, points, weight_attr, decimals=2):
        """Compute and return the weighted mean center (average x and y coordinates weighted by an attribute) of a PointDataset ('points') as a tuple.
        
        Parameters:
        
        points (object): the dataset to analyze. The input dataset will not be altered by the operation.
        
        weight_attr (string): the point attribute to weigh each point's contribution to the mean center by. Must be a numerical attribute present in each point.
        
        decimals (integer): the number of decimal places coordinates will be rounded to. Default is 2."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if len(points) == 0:
            raise ValueError("Cannot compute weighted mean center of an empty dataset.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")

        try:
            weights = [getattr(p, weight_attr) if hasattr(p, weight_attr) else p.attrs[points.headers.index(weight_attr)]
                    for p in points]
        except Exception:
            raise ValueError(f"Attribute '{weight_attr}' not found in points.")

        weight_sum = sum(weights)
        if weight_sum == 0:
            raise ValueError("Sum of weights must be greater than zero.")

        x_sum = sum(p.x * w for p, w in zip(points, weights))
        y_sum = sum(p.y * w for p, w in zip(points, weights))

        return round(x_sum / weight_sum, decimals), round(y_sum / weight_sum, decimals)

    def standard_distance(self, points, decimals=2):
        """Compute and return the standard distance (the average distance of each point from the mean center) of a PointDataset ('points').
        
        Parameters:
        
        points (object): the dataset to analyze. The input dataset will not be altered by the operation.
        
        decimals (integer): the number of decimal places to round to. Default is 2."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if len(points) == 0:
            raise ValueError("Cannot compute standard distance of an empty dataset.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")

        lp = len(points)
        x_mean = sum(p.x for p in points) / lp
        y_mean = sum(p.y for p in points) / lp

        sum_sqr_x = sum((p.x - x_mean) ** 2 for p in points)
        sum_sqr_y = sum((p.y - y_mean) ** 2 for p in points)

        sd = math.sqrt((sum_sqr_x / lp) + (sum_sqr_y / lp))
        return round(sd, decimals)

    def weighted_standard_distance(self, points, weight_attr, decimals=2):
        """Compute and return the weighted standard distance (the average distance of each point from the weighted mean center) of a PointDataset ('points').
        
        Parameters:
        
        points (object): the dataset to analyze. The input dataset will not be altered by the operation.
        
        weight_attr (string): the point attribute to weigh each point's contribution to the mean center by. Must be a numerical attribute present in each point.
        
        decimals (integer): the number of decimal places to round to. Default is 2."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if len(points) == 0:
            raise ValueError("Cannot compute weighted standard distance of an empty dataset.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")

        # Extract weights
        try:
            weights = [getattr(p, weight_attr) if hasattr(p, weight_attr) else p.attrs[points.headers.index(weight_attr)]
                    for p in points]
        except Exception:
            raise ValueError(f"Attribute '{weight_attr}' not found in points.")

        weight_sum = sum(weights)
        if weight_sum == 0:
            raise ValueError("Sum of weights must be greater than zero.")

        # Weighted mean center
        x_mean = sum(p.x * w for p, w in zip(points, weights)) / weight_sum
        y_mean = sum(p.y * w for p, w in zip(points, weights)) / weight_sum

        # Weighted squared deviations
        sum_sqr_x = sum(w * (p.x - x_mean) ** 2 for p, w in zip(points, weights))
        sum_sqr_y = sum(w * (p.y - y_mean) ** 2 for p, w in zip(points, weights))

        sd = math.sqrt((sum_sqr_x + sum_sqr_y) / weight_sum)
        return round(sd, decimals)

    def point_density(self, points, decimals=2):
        """Calculate the density of points (points per unit area) within the extent of a PointDataset ('points').
        
        Parameters:
        
        points (object): the PointDataset the operation will be performed on. The input dataset will not be altered by the operation.
        
        decimals (integer): the number of decimal places the density will be rounded to. Default is 2."""
        if not isinstance(points, _PointDataset):
            raise TypeError("'points' must be a PointDataset object.")
        if len(points) == 0:
            raise ValueError("Cannot calculate density of an empty dataset.")
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("'decimals' must be a non-negative integer.")

        area = self.extent_area(points, decimals=decimals)
        density = len(points) / area
        return round(density, decimals)
     
class _Text():
    """A set of functions for performing operations on lists of string (text) data."""
    # ============================================================
    # PRIVATE METHODS AND CLASS ATTRIBUTES AND FUNCTION DISPLAY
    # ============================================================
    def __init__(self):
        self._total_funcs = len([f for f in dir(__class__) if not f.startswith('__')])
    
    def __repr__(self):
        """If the toolbox is printed, display a message."""
        return ".text: a set of functions for performing operations on lists of string (text) data."

    def display_toolbox_functions(self):
        """Display a list of all available functions within this toolbox."""
        print(f"Number of {__class__.__name__[1:]} functions: {len([f for f in dir(__class__) if not f.startswith('__')])}")
        for f in [f for f in dir(__class__) if not f.startswith("__")]:
            print(f)

    # ============================================================
    # CASE AND CHARACTER TRANSFORMATION
    # ============================================================

    def case(self, vals, mode="upper"):
        """Convert text values from the input list ('vals') to the specified case and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        mode (string): the case conversion mode. Possible values are "upper" to convert to upper case, and "lower" to convert to lower case."""

        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, str) for v in vals):
            raise ValueError("All values in 'vals' must be strings.")

        if mode == "upper":
            return [v.upper() for v in vals]
        elif mode == "lower":
            return [v.lower() for v in vals]
        else:
            raise ValueError("'mode' must be either 'upper' or 'lower'.")
    
    def capitalize(self,vals):
        """Convert the first character of each text value from the input list ('vals ') to upper case and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, str) for v in vals):
            raise ValueError("All values in 'vals' must be strings.")
        return [val.capitalize() for val in vals]

    def reverse(self, vals):
        """Reverse each text value in the input list ('vals') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, str) for v in vals):
            raise ValueError("All values in 'vals' must be strings.")
        return [v[::-1] for v in vals]

    def truncate(self, vals, length, ellipsis=False):
        """Truncate each text value in the input list ('vals') to the specified maximum length ('length'), and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        length (integer): the maximum length for each text value. Must be non-negative.

        ellipsis (Boolean): whether to append '...' if truncation occurs. Defaults to False."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, str) for v in vals):
            raise ValueError("All values in 'vals' must be strings.")

        if not isinstance(length, int) or length < 0:
            raise ValueError("'length' must be a non-negative integer.")
        if not isinstance(ellipsis, bool):
            raise ValueError("'ellipsis' must be a boolean.")
        output = []
        for v in vals:
            if len(v) <= length:
                output.append(v)
            else:
                if ellipsis and length > 3:
                    output.append(v[:length-3] + "...")
                else:
                    output.append(v[:length])
        return output

    # ============================================================
    # CHARACTER AND SUBSTRING ACCESS
    # ============================================================

    def char_at_index(self, vals, index):
        """Find the character at the specified index ('index') of each text value in the input list ('vals') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        index (integer): the index of the character to return from each text value of the input list. If 'index' is not within the valid range for a text value, a value of None will be returned for that element."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, str) for v in vals):
            raise ValueError("All values in 'vals' must be strings.")
        if not isinstance(index, int):
            raise ValueError("'index' must be an integer.")
        return [v[index] if -len(v) <= index < len(v) else None for v in vals]

    def index_of_char(self, vals, char):
        """Find the index of the first occurrence of a specified character ('char') in each text value of the input list ('vals') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        char (string): the case-sensitive character to search for. Must be a non-empty string of length 1. Only the index of the first occurrence will be returned. If the character is not found, a value of None will be returned for that element."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, str) for v in vals):
            raise ValueError("All values in 'vals' must be strings.")
        if not isinstance(char, str):
            raise ValueError("'char' must be a string.")
        if len(char) != 1:
            raise ValueError("'char' must be a single character.")
        return [v.index(char) if char in v else None for v in vals]

    def extract_substring(self, vals, start, end=None):
        """Extract a substring from each text value in the input list ('vals') using the specified start and end indices, and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        start (integer): the starting index (inclusive) of the substring.

        end (integer): the ending index (exclusive) of the substring. If None, the substring extends to the end of the text value."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, str) for v in vals):
            raise ValueError("All values in 'vals' must be strings.")
        if not isinstance(start, int):
            raise ValueError("'start' must be an integer.")
        if end is not None and not isinstance(end, int):
            raise ValueError("'end' must be an integer or None.")
        output = []
        for v in vals:
            try:
                substring = v[start:end]
                output.append(substring if substring != "" else None)
            except Exception:
                output.append(None)
        return output

    def match_position(self, vals, substr, mode="start"):
        """Check whether each text value in the input list ('vals') starts or ends with the specified substring ('substr'), depending on the chosen mode, and return a new list of booleans.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        substr (string): the case-sensitive substring to check for. Must be a non-empty string.

        mode (string): the position to check. Possible values are "start" to check if each text value starts with 'substr', and "end" to check if each text value ends with 'substr'."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, str) for v in vals):
            raise ValueError("All values in 'vals' must be strings.")
        if not isinstance(substr, str):
            raise ValueError("'substr' must be a string.")
        if substr == "":
            raise ValueError("'substr' must be a non-empty string.")
        if mode not in {"start", "end"}:
            raise ValueError("'mode' must be either 'start' or 'end'.")
        if mode == "start":
            return [v.startswith(substr) for v in vals]
        else:  # mode == "end"
            return [v.endswith(substr) for v in vals]

    def find_all(self, vals, substr):
        """Find all indices where the specified substring ('substr') occurs in each text value of the input list ('vals'), and return a nested list of indices.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        substr (string): the case-sensitive substring to search for. Must be a non-empty string."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, str) for v in vals):
            raise ValueError("All values in 'vals' must be strings.")
        if not isinstance(substr, str):
            raise ValueError("'substr' must be a string.")
        if substr == "":
            raise ValueError("'substr' must be a non-empty string.")
        output = []
        for v in vals:
            indices = []
            start = 0
            while True:
                idx = v.find(substr, start)
                if idx == -1:
                    break
                indices.append(idx)
                start = idx + 1
            output.append(indices)
        return output 

    def length(self,vals):
        """Calculate the length of each text value of the input list ('vals') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, str) for v in vals):
            raise ValueError("All values in 'vals' must be strings.")
        return [len(val) for val in vals]

    # ============================================================
    # SEARCHING AND MATCHING
    # ============================================================

    def count_str(self,vals,substr):
        """Count the frequency of occurrence of a specified substring or character ('substr') in each text value of an input list ('vals') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        substr (string): the case-sensitive substring or character to determine the count of."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, str) for v in vals):
            raise ValueError("All values in 'vals' must be strings.")
        if not isinstance(substr, str):
            raise ValueError("'substr' must be a string.")
        if substr == "":
            raise ValueError("'substr' must be a non-empty string.")
        return [val.count(substr) for val in vals]

    def contains(self, vals, substr):
        """Check whether each text value in the input list ('vals') contains the specified substring ('substr') and return a new list of booleans.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        substr (string): the case-sensitive substring to search for. Must be a non-empty string."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, str) for v in vals):
            raise ValueError("All values in 'vals' must be strings.")
        if not isinstance(substr, str):
            raise ValueError("'substr' must be a string.")
        if substr == "":
            raise ValueError("'substr' must be a non-empty string.")
        return [substr in v for v in vals]

    def is_type(self, vals, mode="alpha"):
        """Check whether each text value in the input list ('vals') matches the specified type ('mode'), and return a new list of booleans.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        mode (string, optional): the type of check to perform. Possible values are "alpha" to check if all characters are alphabetic, "digit" to check if all characters are digits, and "alnum" to check if all characters are alphanumeric."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, str) for v in vals):
            raise ValueError("All values in 'vals' must be strings.")
        if mode not in {"alpha", "digit", "alnum"}:
            raise ValueError("'mode' must be 'alpha', 'digit', or 'alnum'.")
        if mode == "alpha":
            return [v.isalpha() for v in vals]
        elif mode == "digit":
            return [v.isdigit() for v in vals]
        else:  # mode == "alnum"
            return [v.isalnum() for v in vals]

    def replace(self,vals,old_str,new_str):
        """Replace a specified subtring or character ('old_str') from each text value of an input list ('vals') with a new specified substring or character ('new_str') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        old_str (string): the old substring or character to be replaced.

        new_str (string): the new substring or character to replace the old."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, str) for v in vals):
            raise ValueError("All values in 'vals' must be strings.")
        if not isinstance(old_str, str):
            raise ValueError("'old_str' must be a string.")
        if old_str == "":
            raise ValueError("'old_str' must be a non-empty string.")
        if not isinstance(new_str, str):
            raise ValueError("'new_str' must be a string.")
        return [val.replace(old_str,new_str) for val in vals]

    # ============================================================
    # CLEANING AND NORMALIZATION
    # ============================================================

    def strip(self, vals):
        """Remove all leading and trailing whitespace from each text value in the input list ('vals') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, str) for v in vals):
            raise ValueError("All values in 'vals' must be strings.")
        return [v.strip() for v in vals]

    def remove_whitespace(self, vals):
        """Remove all whitespace characters from each text value in the input list ('vals') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, str) for v in vals):
            raise ValueError("All values in 'vals' must be strings.")
        return ["".join(v.split()) for v in vals]

    def pad(self, vals, length, char=" ", mode="right"):
        """Pad each text value in the input list ('vals') with the specified character ('char') until it reaches the target length ('length'), and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        length (integer): the target length for each text value.

        char (string): the padding character. Defaults to a space (" "). Must be a single character.

        mode (string, optional): the padding direction. Possible values are "right" to pad on the right side (default), "left" to pad on the left side, and "both" to pad evenly on both sides (extra char goes to the right if uneven)."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, str) for v in vals):
            raise ValueError("All values in 'vals' must be strings.")
        if not isinstance(length, int) or length < 0:
            raise ValueError("'length' must be a non-negative integer.")
        if not isinstance(char, str) or len(char) != 1:
            raise ValueError("'char' must be a single character string.")
        if mode not in {"right", "left", "both"}:
            raise ValueError("'mode' must be 'right', 'left', or 'both'.")
        output = []
        for v in vals:
            if len(v) >= length:
                output.append(v)
            elif mode == "right":
                output.append(v.ljust(length, char))
            elif mode == "left":
                output.append(v.rjust(length, char))
            else:  # mode == "both"
                total_pad = length - len(v)
                left_pad = total_pad // 2
                right_pad = total_pad - left_pad
                output.append(char * left_pad + v + char * right_pad)
        return output

    def trim(self, vals, trim_chars, mode="both"):
        """Trim the specified number of characters ('trim_chars') from text values in the input list ('vals') according to the chosen mode, and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        trim_chars (integer): the number of characters to remove.

        mode (string): the trimming mode. Possible values are "start" to trim characters from the beginning, "end" to trim characters from the end, and "both" to trim characters from both beginning and end."""

        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, str) for v in vals):
            raise ValueError("All values in 'vals' must be strings.")

        if not isinstance(trim_chars, int) or trim_chars < 0:
            raise ValueError("'trim_chars' must be a non-negative integer.")

        if mode not in {"start", "end", "both"}:
            raise ValueError("'mode' must be one of: 'start', 'end', 'both'.")

        if mode == "start":
            return [v[trim_chars:] for v in vals]
        elif mode == "end":
            return [v[:-trim_chars] if trim_chars else v for v in vals]
        else:  # both
            return [v[trim_chars:-trim_chars] if trim_chars else v for v in vals]

    # ============================================================
    # STRUCTURAL OPERATIONS
    # ============================================================

    def split(self, vals, sep=None, maxsplit=-1):
        """Split each text value in the input list ('vals') into a list of substrings, using the specified separator ('sep'), and return a nested list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        sep (string): the delimiter to split on. Defaults to None, which splits on any whitespace.

        maxsplit (integer): the maximum number of splits to perform. Defaults to -1 (no limit)."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, str) for v in vals):
            raise ValueError("All values in 'vals' must be strings.")
        if sep is not None and not isinstance(sep, str):
            raise ValueError("'sep' must be a string or None.")
        if not isinstance(maxsplit, int) or maxsplit < -1:
            raise ValueError("'maxsplit' must be an integer greater than or equal to -1.")

        return [v.split(sep, maxsplit) for v in vals]

    def concatenate(self,vals,list_or_string):
        """Concatenate each text value of the input list ('vals') with either a specified string ('list_or_string') or the correspondong values of a second list ('list_or_string') and return a new list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        list_or_string (string, list): the list of corresponding text values or specified string to concatenate with the input list."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, str) for v in vals):
            raise ValueError("All values in 'vals' must be strings.")
        if isinstance(list_or_string,str):
            return [val + list_or_string for val in vals]
        elif isinstance(list_or_string,list):
            return [val + list_or_string[vals.index(val)] for val in vals]

    def join_strings(self, vals, sep=""):
        """Join all text values in the input list ('vals') into a single string, with an optional separator ('sep') and return a single string.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        sep (string): the separator to insert between each text value. Defaults to an empty string."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, str) for v in vals):
            raise ValueError("All values in 'vals' must be strings.")
        if not isinstance(sep, str):
            raise ValueError("'sep' must be a string.")

        return sep.join(vals)
    
    def to_list(self, vals):
        """Convert each string in the input list ('vals') into a list of characters and return a nested list.

        Parameters:

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if not vals:
            raise ValueError("'vals' must be a non-empty list.")
        if not all(isinstance(v, str) for v in vals):
            raise ValueError("All values in 'vals' must be strings.")
        return [list(v) for v in vals]
   
class _Date():
    """A set of functions for performing operations on lists of dates or lists of numerical values using a corresponding list of dates."""
    # ============================================================
    # PRIVATE METHODS AND CLASS ATTRIBUTES AND FUNCTION DISPLAY
    # ============================================================
    def __init__(self):
        self._total_funcs = len([f for f in dir(__class__) if not f.startswith('__')])
    
    def __repr__(self):
        """If the toolbox is printed, display a message."""
        return ".date: a set of functions for performing operations on lists of dates or lists of numerical values using a corresponding list of dates."

    def display_toolbox_functions(self):
        """Display a list of all available functions within this toolbox."""
        print(f"Number of {__class__.__name__[1:]} functions: {len([f for f in dir(__class__) if not f.startswith('__')])}")
        for f in [f for f in dir(__class__) if not f.startswith("__")]:
            print(f)

    # ============================================================
    # VALIDATION AND NORMALIZATION
    # ============================================================

    def date_to_iso(self,dates,formats):
        """Convert a specified column ('dates') to ISO format (yyyy-mm-dd). Dates that do not match the specified formats ('formats') will be returned unaltered.
        
        Parameters:
        
        dates (list): the list of dates the operation will be performed on. The input list will not be altered by the operation.
        
        formats (string, list): a string or list of strings of date formats to attempt to convert to ISO format."""

        # Normalize formats into a list
        if isinstance(formats, str):
            formats = [formats]
        elif not isinstance(formats, list):
            raise TypeError("formats must be a string or a list of strings.")
        if not all(isinstance(d,str) for d in dates):
            raise ValueError("All dates must be strings.")
        
        strptime = datetime.datetime.strptime
        iso_fmt = "%Y-%m-%d"

        # Prepare list to store converted ISO dates
        output = []

        for val in dates:

            # Default: preserve original value
            iso_val = val

            if isinstance(val, str):
                for fmt in formats:
                    try:
                        dt = strptime(val, fmt)
                        iso_val = dt.strftime(iso_fmt)
                        break
                    except Exception:
                        continue

            output.append(iso_val)

        return output

    def validate_dates(self, dates, raise_errors=True):
        """Validate that all items in a list of dates ('dates') are properly formatted ISO-style date strings. Can return either a Boolean or raise errors.

        Parameters:

        dates (list): the list of dates to validate. Each date must be in a 'yyyy-mm-dd' format. Any character may be used to separate date parts.

        raise_errors (bool): flag to indicate whether to raise detailed ValueError messages when invalid dates are encountered (default True). If False, the function will return True/False instead."""
        if not isinstance(dates, list):
            if raise_errors:
                raise ValueError("'dates' must be a list.")
            return False
        if not all(isinstance(d, str) for d in dates):
            if raise_errors:
                raise ValueError("All dates must be strings.")
            return False

        for d in dates:
            sep = next((ch for ch in d if not ch.isdigit()), None)
            if sep is None:
                if raise_errors:
                    raise ValueError(f"Invalid date format '{d}'. Expected 'yyyy{sep}mm{sep}dd'.")
                return False

            parts = d.split(sep)
            if len(parts) != 3:
                if raise_errors:
                    raise ValueError(f"Invalid date format '{d}'. Expected three parts separated by '{sep}'.")
                return False

            year, month, day = parts
            if not (year.isdigit() and month.isdigit() and day.isdigit()):
                if raise_errors:
                    raise ValueError(f"Invalid date components in '{d}'. Must be numeric.")
                return False
            if len(year) != 4 or len(month) != 2 or len(day) != 2:
                if raise_errors:
                    raise ValueError(f"Invalid date format '{d}'. Expected 'yyyy{sep}mm{sep}dd'.")
                return False
            if not (1 <= int(month) <= 12):
                if raise_errors:
                    raise ValueError(f"Invalid month '{month}' in '{d}'. Must be between 01 and 12.")
                return False
            if not (1 <= int(day) <= 31):
                if raise_errors:
                    raise ValueError(f"Invalid day '{day}' in '{d}'. Must be between 01 and 31.")
                return False

        return True

    def normalize_dates(self, dates, sep="-"):
        """Normalize separators in a list of dates ('dates') to a consistent ISO-style format and return a list of dates.

        Parameters:

        dates (list): the list of dates to normalize. Any character may be used to separate date parts. If a date contains no separators but is in 'yyyymmdd' format, separators will be added automatically.

        sep (string): the separator to use in the normalized output. Default is '-'."""
        if not isinstance(dates, list):
            raise ValueError("'dates' must be a list.")
        if not all(isinstance(d, str) for d in dates):
            raise ValueError("All dates must be strings.")

        output = []
        for d in dates:
            # detect separator by finding the first non-digit character
            sep_in = next((ch for ch in d if not ch.isdigit()), None)

            if sep_in is None:
                # handle raw digit format 'yyyymmdd'
                if len(d) == 8 and d.isdigit():
                    year, month, day = d[:4], d[4:6], d[6:]
                    output.append(f"{year}{sep}{month}{sep}{day}")
                else:
                    raise ValueError(f"Invalid date format '{d}'. Expected 'yyyy{sep}mm{sep}dd'.")
            else:
                parts = d.split(sep_in)
                if len(parts) != 3:
                    raise ValueError(f"Invalid date format '{d}'. Expected three parts separated by '{sep_in}'.")
                year, month, day = parts
                output.append(f"{year}{sep}{month}{sep}{day}")

        return output

    def check_date_continuity(self, dates):
        """Check whether a list of dates ('dates') forms a continuous, gap-free daily sequence from the earliest date to the latest date and return True or False.

        Parameters:

        dates (list): the list of dates the operation will be performed on. The input list will not be altered by the operation. Any character may be used as a separator for string dates. The input list will not be altered by the operation."""
        # Basic type validation
        if not isinstance(dates, list) or not dates:
            raise ValueError("'dates' must be a non-empty list.")
        if not all(isinstance(d,str) for d in dates):
            raise ValueError("All dates must be strings.")

        # sort dates and get range
        sorted_dates = sorted(dates)
        start, end = self.get_date_range(sorted_dates)

        # generate_day_range is non-inclusive of end, so add one day manually
        end_dt = datetime.datetime.strptime(end, "%Y-%m-%d") + datetime.timedelta(days=1)
        end_plus = end_dt.strftime("%Y-%m-%d")

        full_range = self.generate_day_range(start, end_plus)

        return len(full_range) == len(sorted_dates)

    def find_missing_dates(self, dates):
        """Return a sorted list of missing dates within the range defined by the earliest and latest dates in the input list ('dates').

        Parameters:

        dates (list): the list of dates the operation will be performed on. All dates must be ISO-style strings ('yyyy-mm-dd'). The input list will not be altered by the operation."""
        # Basic type validation
        if not isinstance(dates, list) or not dates:
            raise ValueError("'dates' must be a non-empty list.")
        if not all(isinstance(d, str) for d in dates):
            raise ValueError("All dates must be strings.")

        # Sort dates and get range
        sorted_dates = sorted(dates)
        start, end = self.get_date_range(sorted_dates)

        # generate_day_range is non-inclusive of end, so add one day manually
        end_dt = datetime.datetime.strptime(end, "%Y-%m-%d") + datetime.timedelta(days=1)
        end_plus = end_dt.strftime("%Y-%m-%d")

        full_range = self.generate_day_range(start, end_plus)

        # Identify missing dates
        date_set = set(sorted_dates)
        missing = [d for d in full_range if d not in date_set]

        return missing

    def insert_missing_date_values(self, dates, vals, placeholder=""):
        """ Insert placeholder values into a list of values ('vals') wherever dates are missing from the list of dates ('dates'). Returns both the new date list and a new value list with placeholders inserted.

        Parameters:

        dates (list): the list of ISO-style dates ('yyyy-mm-dd') the operation will be performed on. The input list will not be altered.

        vals (list): the list of values corresponding to each date. Must be the same length as 'dates'.

        placeholder: the value to insert for missing dates. Default is ""."""
        # Basic type validation
        if not isinstance(dates, list) or not dates:
            raise ValueError("'dates' must be a non-empty list.")
        if not isinstance(vals, list):
            raise ValueError("'vals' must be a list.")
        if len(dates) != len(vals):
            raise ValueError("'dates' and 'vals' must be the same length.")
        if not all(isinstance(d, str) for d in dates):
            raise ValueError("All dates must be strings.")

        # Sort dates and values together 
        paired = sorted(zip(dates, vals), key=lambda x: x[0])
        sorted_dates = [d for d, _ in paired]
        sorted_values = [v for _, v in paired]

        # Determine full expected date range
        start, end = self.get_date_range(sorted_dates)

        # generate_day_range is non-inclusive of end, so add one day manually
        end_dt = datetime.datetime.strptime(end, "%Y-%m-%d") + datetime.timedelta(days=1)
        end_plus = end_dt.strftime("%Y-%m-%d")

        full_range = self.generate_day_range(start, end_plus)

        # Build new aligned lists
        date_set = set(sorted_dates)
        value_map = dict(zip(sorted_dates, sorted_values))

        new_dates = []
        new_values = []

        for d in full_range:
            new_dates.append(d)
            if d in date_set:
                new_values.append(value_map[d])
            else:
                new_values.append(placeholder)

        return new_dates, new_values

    def is_leap_year(self,year):
        """Determine if the specified year is a leap year. Returns True or False.
        
        Parameters:
        
        year (integer, string): the year the operation will be performed on."""
        if not isinstance(year, int) and not isinstance(year, str):
            raise ValueError("Year value must be an integer or a string.")
        if isinstance(year,str):
            year = int(year)
        
        return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)

    # ============================================================
    # DATE RANGE GENERATION
    # ============================================================
    
    def generate_day_range(self,start,end):
        """Return a list of daily dates in a yyyy-mm-dd format within the specified range, with years starting on January 1 (01-01).

        Parameters:

        start (string): the starting date in the date range. Must be in a "yyyy-mm-dd" format. Any character may be used to separate date parts, but separators must be the same between start and end dates.
        
        end (string): the ending date in the date range. This date is non-inclusive (eg. if an ending date of 2020-01-01 is specified, the last day in the output date range will be 2019-12-31). Must be in a "yyyy-mm-dd" format. A hyphen "-", forward slash "/", or space " " may be used to separate date parts."""
        # detect separator by finding the first non-digit character
        sep_start = next((ch for ch in start if not ch.isdigit()), None)
        sep_end = next((ch for ch in end if not ch.isdigit()), None)

        if sep_start is None or sep_end is None:
            raise ValueError("Dates must contain a separator between year, month, and day.")
        if sep_start != sep_end:
            raise ValueError("Start and end dates must use the same separator.")
        if not isinstance(start, str) or not isinstance(end, str):
            raise ValueError("Start and end dates must be strings in 'yyyy-mm-dd' format.")

        d_sep = sep_start
        start_dt = datetime.datetime.strptime(start, f"%Y{d_sep}%m{d_sep}%d")
        end_dt = datetime.datetime.strptime(end, f"%Y{d_sep}%m{d_sep}%d")

        date_generated = [start_dt + datetime.timedelta(days=x) for x in range((end_dt - start_dt).days)]
        return [d.strftime(f"%Y{d_sep}%m{d_sep}%d") for d in date_generated]
    
    def generate_day_range_julian(self,start,end):
        """Return a list of daily dates in a Julian format within the specified range, with years starting on January 1 (001).

        Parameters:

        start (string): the starting date in the date range. Must be in a "yyyy-mm-dd" format. Any character may be used to separate date parts.
        
        end (string): the ending date in the date range. This date is non-inclusive (eg. if an ending date of 2020-01-01 is specified, the last day in the output date range will be 2019-12-31 (2019365)). Must be in a "yyyy-mm-dd" format. A hyphen "-", forward slash "/", or space " " may be used to separate date parts."""
        # detect separator by finding the first non-digit character
        sep_start = next((ch for ch in start if not ch.isdigit()), None)
        sep_end = next((ch for ch in end if not ch.isdigit()), None)
        if sep_start is None or sep_end is None or sep_start != sep_end:
            raise ValueError("Start and end dates must use the same separator.")
        if not isinstance(start, str) or not isinstance(end, str):
            raise ValueError("Start and end dates must be strings in 'yyyy-mm-dd' format.")

        d_sep = sep_start
        start_dt = datetime.datetime.strptime(start, f"%Y{d_sep}%m{d_sep}%d")
        end_dt = datetime.datetime.strptime(end, f"%Y{d_sep}%m{d_sep}%d")

        date_generated = [start_dt + datetime.timedelta(days=x) for x in range((end_dt - start_dt).days)]
        return [f"{d.year}{d.timetuple().tm_yday:03d}" for d in date_generated]

    def generate_day_range_doy(self,start,end):
        """Return a list of daily dates in a day-of-year format within the specified range.

        Parameters:

        start (string): the starting date in the date range. Must be in a "yyyy-mm-dd" format. Any character may be used to separate date parts, but separators must be the same between start and end dates.
        
        end (string): the ending date in the date range. This date is non-inclusive (eg. if an ending date of 2020-01-01 is specified, the last day in the output date range will be 2019-12-31). Must be in a "yyyy-mm-dd" format. A hyphen "-", forward slash "/", or space " " may be used to separate date parts."""
        # detect separator by finding the first non-digit character
        sep_start = next((ch for ch in start if not ch.isdigit()), None)
        sep_end = next((ch for ch in end if not ch.isdigit()), None)

        if sep_start is None or sep_end is None:
            raise ValueError("Dates must contain a separator between year, month, and day.")
        if sep_start != sep_end:
            raise ValueError("Start and end dates must use the same separator.")
        if not isinstance(start, str) or not isinstance(end, str):
            raise ValueError("Start and end dates must be strings in 'yyyy-mm-dd' format.")

        d_sep = sep_start
        start_dt = datetime.datetime.strptime(start, f"%Y{d_sep}%m{d_sep}%d")
        end_dt = datetime.datetime.strptime(end, f"%Y{d_sep}%m{d_sep}%d")

        date_generated = [start_dt + datetime.timedelta(days=x) for x in range((end_dt - start_dt).days)]
        return [d.timetuple().tm_yday for d in date_generated]

    def generate_day_range_alt(self,start,end,month):
        """Return a list of daily dates in a yyyy-mm-dd format within the specified range, where the years start at an alternative month (i.e. not January (01)). This function could be used to generate dates by hydrological year (October to September), for example.

        Parameters:

        start (string): the starting date in the date range. Must be in a "yyyy-mm-dd" format. Any character may be used to separate date parts.
        
        end (string): the ending date in the date range. This date is non-inclusive (eg. if an ending date of 2020-01-01 is specified, the last day in the output date range will be 2019-12-31). Must be in a "yyyy-mm-dd" format. A hyphen "-", forward slash "/", or space " " may be used to separate date parts.
        
        month (integer): the integer corresponding to the month each year will begin at. Must be greater than 1 (January) and less than or equal to 12 (December)."""
        if not isinstance(start, str) or not isinstance(end, str):
            raise ValueError("Start and end dates must be strings in 'yyyy-mm-dd' format.")
        if not isinstance(month, int) or month < 2 or month > 12:
            raise ValueError("Month must be an integer between 2 and 12.")

        years = self.generate_day_range(start, end)
        output = []
        for y in years:
            if int(y[5:7]) >= month:
                y = str(int(y[:4]) + 1) + y[4:]
            output.append(y)
        return output

    def generate_month_range(self,start = 1,length = 12,form = "numeric"):
        """Return a list of months of the specified length ('length') and the specified format ('format').

        Parameters:

        start (integer): the integer corresponding to the starting month in the month range. Should be greater than or equal to 1 (January) and less than or equal to 12 (December)

        length (integer): the number of months in the output range.

        form (string): the format of months in the output list. May be specified as "numeric" for integer format (1,2,3...), "padded numeric" for integers as zero-padded strings ("01","02","03"...), "full text" for full month names ("January","February","March"...) or "short text" for shortened month names ("Jan","Feb","Mar"...)."""
        if not isinstance(start, int) or start < 1 or start > 12:
            raise ValueError("Start must be an integer between 1 and 12.")
        if not isinstance(length, int) or length < 0:
            raise ValueError("Length must be a non-negative integer.")
        valid_forms = {"numeric", "padded numeric", "full text", "short text"}
        if not isinstance(form, str) or form not in valid_forms:
            raise ValueError("Form must be one of: 'numeric', 'padded numeric', 'full text', 'short text'.")

        # Base sequences
        numeric = [i for i in range(1, 13)]
        padded_numeric = [str(m).zfill(2) for m in numeric]
        full_text = ["January","February","March","April","May","June","July","August","September","October","November","December"]
        short_text = [m[:3] for m in full_text]

        # Choose formatter
        formats = {
            "numeric": numeric,
            "padded numeric": padded_numeric,
            "full text": full_text,
            "short text": short_text,
        }
        seq = formats[form]

        # Build output with wrapping
        output = []
        idx = start - 1
        for _ in range(length):
            output.append(seq[idx])
            idx += 1
            if idx >= 12:
                idx = 0
        return output

    def generate_year_range(self,start,end):
        """Return a list of yearly dates in a yyyy format within the specified range.

        Parameters:

        start (integer, string): the starting date in the date range. May be specified as an integer in yyyy format or an ISO format date string.
        
        end (integer ,string): the ending date in the date range. This date is non-inclusive (eg. if an ending date of 2020 is specified, the last year in the output date range will be 2019). May be specified as an integer in yyyy format or an ISO format date string."""
        if isinstance(start,str):
            start = int(start[:4])
        if isinstance(end,str):
            end = int(end[:4])
        if start >= end:
            raise ValueError("Start year must be less than end year.")
        return list(range(start, end))
           
    def generate_leap_years(self,start,end):
        """Return a list of leap years in a yyyy format within the specified range.

        Parameters:

        start (integer): the starting date in the date range. Must be in a yyyy format.
        
        end (integer): the ending date in the date range. This date is non-inclusive (eg. if an ending date of 2020 is specified, the last year in the output date range will be 2019). Must be in a yyyy format."""
        if not isinstance(start, int) or not isinstance(end, int):
            raise ValueError("Start and end must be integers representing years (e.g., 2020).")
        if start >= end:
            raise ValueError("Start year must be less than end year.")

        return [d for d in range(start,end) if d % 4 == 0 and (d % 100 != 0 or d % 400 == 0)]

    def extend_dates(self,dates,num_days):
        """Extend the input list of dates ('dates') by a specified number of days ('num_days') and return a new list of dates.
        
        Parameters:
        
        dates (list): the list of dates the operation will be performed on. Must be in a "yyyy-mm-dd" format. Any character may be used to separate date parts. The input list will not be altered by the operation.
        
        num_days (integer): the number of days to extend the list of dates by."""
        # Validate inputs
        if not isinstance(dates, list) or not dates:
            raise ValueError("Dates must be a non-empty list of date strings.")
        if not all(isinstance(d, str) for d in dates):
            raise ValueError("All dates must be strings in 'yyyy-mm-dd' format.")
        if not isinstance(num_days, int) or num_days < 0:
            raise ValueError("num_days must be a non-negative integer.")

        # Detect separator from first date
        sep_first = next((ch for ch in dates[0] if not ch.isdigit()), None)
        if sep_first is None:
            raise ValueError("Dates must contain a separator between year, month, and day.")

        # Ensure all dates use the same separator
        for d in dates:
            sep_d = next((ch for ch in d if not ch.isdigit()), None)
            if sep_d != sep_first:
                raise ValueError("All dates must use the same separator.")

        d_sep = sep_first
        date_copy = dates[:]

        # Parse last date in list
        last_dt = datetime.datetime.strptime(date_copy[-1], f"%Y{d_sep}%m{d_sep}%d")

        # Extend by num_days
        for d in range(1, num_days + 1):
            new_dt = last_dt + datetime.timedelta(days=d)
            date_copy.append(new_dt.strftime(f"%Y{d_sep}%m{d_sep}%d"))

        return date_copy

    # ============================================================
    # DATE RANGE OPERATIONS
    # ============================================================    
    
    def get_date_range(self,vals):
        """Return a tuple of (start date, end date) from the input list ('vals') of dates.
        
        Parameters:
        
        vals (list): the list of dates the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(vals, list):
            raise ValueError("Input must be a list of dates.")
        if len(vals) < 2:
            raise ValueError("Date list must contain at least two elements.")

        return min(vals),max(vals)

    def get_years(self,vals,dtype="integer"):
        """Separate the years from each date in an input list of dates ('vals') and return a new list of years.
        
        Parameters:
        
        vals (list): the list of dates the operation will be performed on. The input list will not be altered by the operation. List of dates must be in a 'yyyy-mm-dd' format. Any character may be used to separate date parts.
        
        dtype (string): the data type of the output years. May be specified as "integer" or "string" for integer or string representations of years, respectively."""
        if not isinstance(vals, list):
            raise ValueError("Input must be a list of date strings.")
        if any(len(d) < 4 or not isinstance(d, str) for d in vals):
            raise ValueError("All dates must be strings in 'yyyy-mm-dd' format")
        
        if dtype == "integer":
            return [int(d[:4]) for d in vals]
        elif dtype == "string":
            return [d[:4] for d in vals]
        else:
            raise ValueError("dtype must be either 'integer' or 'string'.")

    def get_months(self,vals,dtype="integer"):
        """Separate the months from each date in an input list of dates ('vals') and return a new list of months.
        
        Parameters:
        
        vals (list): the list of dates the operation will be performed on. The input list will not be altered by the operation. List of dates must be in a 'yyyy-mm-dd' format. Any character may be used to separate date parts.
        
        dtype (string): the data type of the output months. May be specified as "integer" or "string" for integer or string representations of months, respectively."""
        if not isinstance(vals, list):
            raise ValueError("Input must be a list of date strings.")
        if any(len(d) < 7 or not isinstance(d, str) for d in vals):
            raise ValueError("All dates must be strings in 'yyyy-mm-dd' format")
        
        if dtype == "integer":
            return [int(d[5:7]) for d in vals]
        elif dtype == "string":
            return [d[5:7] for d in vals]
        else:
            raise ValueError("dtype must be either 'integer' or 'string'.")         
    
    def get_days(self,vals,dtype="integer"):
        """Separate the days from each date in an input list of dates ('vals') and return a new list of days.
        
        Parameters:
        
        vals (list): the list of dates the operation will be performed on. The input list will not be altered by the operation. List of dates must be in a 'yyyy-mm-dd' format. Any character may be used to separate date parts.
        
        dtype (string): the data type of the output days. May be specified as "integer" or "string" for integer or string representations of days, respectively."""
        if not isinstance(vals, list):
            raise ValueError("Input must be a list of date strings.")
        if any(len(d) < 10 or not isinstance(d, str) for d in vals):
            raise ValueError("All dates must be strings in 'yyyy-mm-dd' format")
        
        if dtype == "integer":
            return [int(d[8:]) for d in vals]
        elif dtype == "string":
            return [d[8:] for d in vals]
        else:
            raise ValueError("dtype must be either 'integer' or 'string'.")       
   
    def date_range_intersection(self,range1,range2):
        """Determine if two input date ranges ('range1','range2') overlap and return the start and end dates of the overlapping range. It is assumed that the date ranges are sorted in order. Returns a tuple of (start, end). If there is no overlap between date ranges, a tuple of (None, None) is returned.
        
        Parameters:
        
        range1 (list): the first list of dates the operation will be performed on. Must be in a "yyyy-mm-dd" format. Any character may be used to separate date parts. The input list will not be altered by the operation.
        
        range1 (list): the second list of dates the operation will be performed on. Must be in a "yyyy-mm-dd" format. Any character may be used to separate date parts. The input list will not be altered by the operation."""
        if not isinstance(range1, list) or not isinstance(range2, list) or not range1 or not range2:
            raise ValueError("Both ranges must be non-empty lists of date strings.")
        if not all(isinstance(d, str) for d in range1 + range2):
            raise ValueError("All dates must be strings in 'yyyy-mm-dd' format.")

        # Detect separators
        sep1 = next((ch for ch in range1[0] if not ch.isdigit()), None)
        sep2 = next((ch for ch in range2[0] if not ch.isdigit()), None)
        if sep1 is None or sep2 is None or sep1 != sep2:
            raise ValueError("Both ranges must use the same separator.")
        d_sep = sep1

        # Extract earliest/latest start/end using raw string comparison
        latest_start = max(range1[0], range2[0])
        earliest_end = min(range1[-1], range2[-1])

        # Check for overlap
        if earliest_end < latest_start:
            return (None, None)
        else:
            return (latest_start, earliest_end)
    
    def date_range_slice(self,dates,vals,start,end):
        """Return a subset of values from an input list ('vals') within a specified date range using a corresponding list of dates ('dates').

        Parameters:

        dates (list): the list of dates used to aggregate the input list of values. List of dates must be in a 'yyyy-mm-dd' format. Any character may be used to separate date parts. List of dates must be the same length as the list of values. It is expected that dates are ordered sequentially.

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        start (string): the starting date in the output range of values. The string "start" may be specified to indicate the start of the date range, otherwise it must be in a "yyyy-mm-dd" format. A hyphen "-", forward slash "/", or space " " may be used to separate date parts.
        
        end (string): the ending date in the output range of values. This date is non-inclusive (eg. if an ending date of 2020-01-01 is specified, the last day in the output date range will be 2019-12-31). The string "end" may be specified to indicate the end of the date range, otherwise it must be in a "yyyy-mm-dd" format. A hyphen "-", forward slash "/", or space " " may be used to separate date parts."""
        if not isinstance(dates, list) or not isinstance(vals, list):
            raise ValueError("Both 'dates' and 'vals' must be lists.")
        if len(dates) != len(vals):
            raise ValueError("'dates' and 'vals' must be the same length.")
        if not all(isinstance(d, str) for d in dates):
            raise ValueError("All dates must be strings in 'yyyy-mm-dd' format.")
        if not isinstance(start, str) or not isinstance(end, str):
            raise ValueError("'start' and 'end' must be strings.")

        # Resolve start index
        if start == "start":
            start_idx = 0
        else:
            try:
                start_idx = dates.index(start)
            except ValueError:
                raise ValueError(f"Start date '{start}' not found in 'dates'.")

        # Resolve end index
        if end == "end":
            end_idx = len(dates)
        else:
            try:
                end_idx = dates.index(end)
            except ValueError:
                raise ValueError(f"End date '{end}' not found in 'dates'.")

        # Slice values
        return vals[start_idx:end_idx]
    
    def threshold_date(self,dates,condition,threshold):
        """Apply a date threshold ('threshold') to a list of dates ('dates') and return a new list.
        
        Parameters:
        
        dates (list): the list of dates used to aggregate the input list of values. List of dates must be in a 'yyyy-mm-dd' format. Any character may be used to separate date parts. List of dates must be the same length as the list of values. It is expected that dates are ordered sequentially.

        condition (string): the condition placed upon each date in the input list. May be specified as "<" for less than, "<=" for less than equal to, ">" for greater than, ">=" for greater than equal to, "==" for equal to, or "!=" for not equal to.
        
        threshold (string): the threshold to evaluate each date in the input list against. Must be in "yyyy-mm-dd" format."""
        if not isinstance(dates, list):
            raise ValueError("'dates' must be a list.")
        if not all(isinstance(d, str) for d in dates):
            raise ValueError("All dates must be strings in 'yyyy-mm-dd' format.")
        if not isinstance(condition, str):
            raise ValueError("'condition' must be a string.")
        if not isinstance(threshold, str):
            raise ValueError("'threshold' must be a string in 'yyyy-mm-dd' format.")

        valid_conditions = {'<', '<=', '>', '>=', '==', '!='}
        if condition not in valid_conditions:
            raise ValueError(f"Unknown condition '{condition}'.")

        # Normalize separators to '-' for consistent lexicographic comparison
        def norm(s):
            return s.replace('/', '-').replace(' ', '-')

        thr = norm(threshold)
        ops = {'<': operator.lt, '<=': operator.le, '>': operator.gt,
               '>=': operator.ge, '==': operator.eq, '!=': operator.ne}

        # Filter using lexicographic comparison of normalized ISO strings
        output = [d for d in dates if ops[condition](norm(d), thr)]
        return output

    def optimal_overlap_range(self,date_ranges, tolerance=0.05):
        """Determine the 'optimal' overlapping date range across many input date ranges ('dates') by evaluating all combinations of start and end dates from the input dates and assigning a score based on number of overlapping date ranges and the total days captured by each combination. Optimal is defined as the start and end date that maximizes the number of input ranges that overlap within the given tolerance of days while favouring longer shared durations. Returns a dictionary of "best_start" (i.e., optimal start date), "best_end" (i.e., optimal end date), "count" (i.e, the number of input date ranges that overlap the optimal range within the tolerance), "percent" (i.e, the percentage of input date ranges that overlap the optimal range within the tolerance), and "days" (i.e., the number of days in the optimal date range). 

        Parameters:

        date_sequences (list): a nested list of date range lists. Individual dates must be in ISO format (i.e., yyyy-mm-dd)
            
        tolerance (float): the tolerance threshold for required overlap as a percent of shared days between date ranges."""
        if not isinstance(date_ranges,list):
            raise ValueError("'date_ranges' must be a list.")
        for r in date_ranges:
            if not isinstance(r, list):
                raise ValueError("Each date range in 'date_ranges' must be a list.")
            if not all(isinstance(d,str) for d in r):
                raise ValueError("Each date range in 'date_ranges' must be a list of ISO strings.")
        if not isinstance(tolerance,float):
            raise ValueError("'tolerance' must be a floating point value")
        if tolerance < 0 or tolerance > 1:
            raise ValueError("'tolerance' must be between 0.0 and 1.0 (inclusive).")
        starts = [datetime.datetime.fromisoformat(min(seq)) for seq in date_ranges]
        ends   = [datetime.datetime.fromisoformat(max(seq)) for seq in date_ranges]

        results = []

        for s_dt in starts:
            for e_dt in ends:
                if s_dt >= e_dt:
                    continue

                target_seconds = (e_dt - s_dt).total_seconds()
                duration_days = (e_dt - s_dt).days

                count = 0

                for seq_start, seq_end in zip(starts, ends):
                    overlap_start = max(s_dt, seq_start)
                    overlap_end   = min(e_dt, seq_end)

                    if overlap_start < overlap_end:
                        overlap_seconds = (overlap_end - overlap_start).total_seconds()
                        if overlap_seconds / target_seconds >= (1 - tolerance):
                            count += 1

                # use an internal score so that it's not just the number of date ranges, but the total days of each
                score = count * duration_days

                results.append({
                    "start": s_dt,
                    "end": e_dt,
                    "count": count,
                    "days": duration_days,
                    "score": score,
                })

        best = max(results, key=lambda x: x["score"])

        percent = (best["count"] / len(date_ranges)) * 100

        return {
            "best_start": best["start"].date().isoformat(), #start date of the optimal date range
            "best_end": best["end"].date().isoformat(), #end date of the optimal date range
            "count": best["count"], #number of input date ranges captured by the optimal range
            "percent": percent, #percentage of input date ranges captured by the optimal range
            "days": best["days"], #number of days in optimal date range
        }

    # ============================================================
    # TEMPORAL AGGREGATION
    # ============================================================
    
    def aggregate_year(self,dates,vals,method):
        """Aggregate the values of an input list ('vals') by year using a corresponding list of dates ('dates') and return a new list.

        Parameters:

        dates (list): the list of dates used to aggregate the input list of values. List of dates must be in a 'yyyy-mm-dd' format. Any character may be used to separate date parts. List of dates must be the same length as the list of values. It is expected that dates are ordered sequentially.

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        method (string): the method used to aggregate the data. Possible aggregation methods are "Count", "Unique Values", "First", "Last", "Sum", "Max", "Min", "Range", "Mean", "Median", "Mode", "IQ1", "IQ3", "IQR", "Standard Deviation", "Standard Error", "Variance", "Skewness", and "Kurtosis" """
        # Validate inputs
        if not isinstance(dates, list) or not isinstance(vals, list):
            raise ValueError("Both 'dates' and 'vals' must be lists.")
        if len(dates) != len(vals):
            raise ValueError("'dates' and 'vals' must be the same length.")
        if not all(isinstance(d, str) for d in dates):
            raise ValueError("All dates must be strings in 'yyyy-mm-dd' format.")
        if not isinstance(method, str):
            raise ValueError("'method' must be a string.")
        valid_stats = {
        "Count", "Unique Values", "First", "Last", "Sum", "Max", "Min", "Range", "Mean",
        "Median", "Mode", "IQ1", "IQ3", "IQR", "Standard Deviation",
        "Standard Error", "Variance", "Coefficient of Variation",
        "Skewness", "Kurtosis"}

        if method not in valid_stats:
            raise ValueError(f"Unknown statistic '{method}'.")
        v = self.group_by_year(dates, vals)
        stats = _Stats()
        output = []
        for d in v:
            output.append(stats.descriptive_stat(d, method))
        return output

    def aggregate_year_single(self,dates,vals,year,method):
        """Aggregate the values of an input list ('vals') for a single specified year ('year') using a corresponding list of dates ('dates') and return the aggregated value.

        Parameters:

        dates (list): the list of dates used to aggregate the input list of values. List of dates must be in a 'yyyy-mm-dd' format. Any character may be used to separate date parts. List of dates must be the same length as the list of values. It is expected that dates are ordered sequentially.

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        year (integer): the single year of data to aggregate. Must be in a 'yyyy' format.
        
        method (string): the method used to aggregate the data. Possible aggregation methods are "Count", "Unique Values", "First", "Last", "Sum", "Max", "Min", "Range", "Mean", "Median", "Mode", "IQ1", "IQ3", "IQR", "Standard Deviation", "Standard Error", "Variance", "Skewness", and "Kurtosis" """
        if not isinstance(dates, list) or not isinstance(vals, list):
            raise ValueError("Both 'dates' and 'vals' must be lists.")
        if len(dates) != len(vals):
            raise ValueError("'dates' and 'vals' must be the same length.")
        if not all(isinstance(d, str) for d in dates):
            raise ValueError("All dates must be strings in 'yyyy-mm-dd' format.")
        if not isinstance(year, int):
            raise ValueError("'year' must be an integer in 'yyyy' format.")
        if not isinstance(method, str):
            raise ValueError("'method' must be a string.")

        valid_stats = {
        "Count", "Unique Values", "First", "Last", "Sum", "Max", "Min", "Range", "Mean",
        "Median", "Mode", "IQ1", "IQ3", "IQR", "Standard Deviation",
        "Standard Error", "Variance", "Coefficient of Variation",
        "Skewness", "Kurtosis"}

        if method not in valid_stats:
            raise ValueError(f"Unknown statistic '{method}'.")
        out = [v for d, v in zip(dates, vals) if int(d[:4]) == year]

        if method == "Count":
            return len(out)
        elif method == "Sum":
            return sum(out)
        else:
            stats = _Stats()
            return stats.descriptive_stat(out, method)

    def aggregate_month(self,dates,vals,method):
        """Aggregate the values of an input list ('vals') by month using a corresponding list of dates ('dates') and return a new list.

        Parameters:

        dates (list): the list of dates used to aggregate the input list of values. List of dates must be in a 'yyyy-mm-dd' format. Any character may be used to separate date parts. List of dates must be the same length as the list of values. It is expected that dates are ordered sequentially.

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        method (string): the method used to aggregate the data. Possible aggregation methods are "Count", "Unique Values", "First", "Last", "Sum", "Max", "Min", "Range", "Mean", "Median", "Mode", "IQ1", "IQ3", "IQR", "Standard Deviation", "Standard Error", "Variance", "Skewness", and "Kurtosis" """
        if not isinstance(dates, list) or not isinstance(vals, list):
            raise ValueError("Both 'dates' and 'vals' must be lists.")
        if len(dates) != len(vals):
            raise ValueError("'dates' and 'vals' must be the same length.")
        if not all(isinstance(d, str) for d in dates):
            raise ValueError("All dates must be strings in 'yyyy-mm-dd' format.")
        if not isinstance(method, str):
            raise ValueError("'method' must be a string.")

        valid_stats = {
        "Count", "Unique Values", "Sum", "Max", "Min", "Range", "Mean",
        "Median", "Mode", "IQ1", "IQ3", "IQR", "Standard Deviation",
        "Standard Error", "Variance", "Coefficient of Variation",
        "Skewness", "Kurtosis"}

        if method not in valid_stats:
            raise ValueError(f"Unknown statistic '{method}'.")
        v = self.group_by_month(dates,vals)
        output = []
        for d in v:
            if method == "First":
                output.append(d[0])
            elif method == "Last":
                output.append(d[-1])
            else:
                output.append(_Stats().descriptive_stat(d,method))
        return output

    def aggregate_season(self,dates,vals,months,method):
        """Aggregate the values of an input list ('vals') by a specified season ('months') using a corresponding list of dates ('dates') and return a new list.

        Parameters:

        dates (list): the list of dates used to aggregate the input list of values. List of dates must be in a 'yyyy-mm-dd' format. Any character may be used to separate date parts. List of dates must be the same length as the list of values. It is expected that dates are ordered sequentially.

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        months (list): a list of integers corresponding to the months within a particular season. For example, if aggregating within the winter season, [1,2,12] would correspond to the months of January, February, and December, respectively. Month values should be given in ascending order.
        
        method (string): the method used to aggregate the data. Possible aggregation methods are "Count", "Unique Values", "First", "Last", "Sum", "Max", "Min", "Range", "Mean", "Median", "Mode", "IQ1", "IQ3", "IQR", "Standard Deviation", "Standard Error", "Variance", "Skewness", and "Kurtosis" """
        if not isinstance(dates, list) or not isinstance(vals, list):
            raise ValueError("Both 'dates' and 'vals' must be lists.")
        if len(dates) != len(vals):
            raise ValueError("'dates' and 'vals' must be the same length.")
        if not all(isinstance(d, str) for d in dates):
            raise ValueError("All dates must be strings in 'yyyy-mm-dd' format.")
        if not isinstance(months, list) or not all(isinstance(m, int) for m in months):
            raise ValueError("'months' must be a list of integers.")
        if not isinstance(method, str):
            raise ValueError("'method' must be a string.")

        valid_stats = {
            "Count", "Unique Values", "First", "Last", "Sum", "Max", "Min", "Range", "Mean",
            "Median", "Mode", "IQ1", "IQ3", "IQR", "Standard Deviation",
            "Standard Error", "Variance", "Coefficient of Variation",
            "Skewness", "Kurtosis"
        }

        if method not in valid_stats:
            raise ValueError(f"Unknown statistic '{method}'.")

        # Grouping
        v = self.group_by_season(dates, vals, months)

        # Aggregation
        stats = _Stats()
        output = []
        for d in v:
            if method == "First":
                output.append(d[0])
            elif method == "Last":
                output.append(d[-1])
            elif method == "Count":
                output.append(len(d))
            elif method == "Sum":
                output.append(sum(d))
            else:
                output.append(stats.descriptive_stat(d, method))
        return output
  
    def aggregate_day(self,dates,vals,method):
        """Aggregate the values of an input list ('vals') by day using a corresponding list of dates ('dates') and return a new list.

        Parameters:

        dates (list): the list of dates used to aggregate the input list of values. List of dates must be in a 'yyyy-mm-dd' format. Any character may be used to separate date parts. List of dates must be the same length as the list of values. It is expected that dates are ordered sequentially.

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        method (string): the method used to aggregate the data. Possible aggregation methods are "Count", "Unique Values", "First", "Last", "Sum", "Max", "Min", "Range", "Mean", "Median", "Mode", "IQ1", "IQ3", "IQR", "Standard Deviation", "Standard Error", "Variance", "Skewness", and "Kurtosis" """
        if not isinstance(dates, list) or not isinstance(vals, list):
            raise ValueError("Both 'dates' and 'vals' must be lists.")
        if len(dates) != len(vals):
            raise ValueError("'dates' and 'vals' must be the same length.")
        if not all(isinstance(d, str) for d in dates):
            raise ValueError("All dates must be strings in 'yyyy-mm-dd' format.")
        if not isinstance(method, str):
            raise ValueError("'method' must be a string.")

        valid_stats = {
            "Count", "Unique Values", "First", "Last", "Sum", "Max", "Min", "Range", "Mean",
            "Median", "Mode", "IQ1", "IQ3", "IQR", "Standard Deviation",
            "Standard Error", "Variance", "Coefficient of Variation",
            "Skewness", "Kurtosis"
        }

        if method not in valid_stats:
            raise ValueError(f"Unknown statistic '{method}'.")

        # Grouping
        v = self.group_by_day(dates, vals)

        # Aggregation
        stats = _Stats()
        output = []
        for d in v:
            if method == "First":
                output.append(d[0])
            elif method == "Last":
                output.append(d[-1])
            elif method == "Count":
                output.append(len(d))
            elif method == "Sum":
                output.append(sum(d))
            else:
                output.append(stats.descriptive_stat(d, method))
        return output

    # ============================================================
    # TEMPORAL GROUPING
    # ============================================================

    def group_by_year(self,dates,vals):
        """Group the values from an input list ('vals') using a corresponding list of dates ('dates') into lists of values within each year in the date range and return a nested list of yearly sublists.

        Parameters:

        dates (list): the list of dates used to aggregate the input list of values. List of dates must be in a 'yyyy-mm-dd' format. Any character may be used to separate date parts. List of dates must be the same length as the list of values. It is expected that dates are ordered sequentially.

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(dates, list) or not isinstance(vals, list):
            raise ValueError("Both 'dates' and 'vals' must be lists.")
        if len(dates) != len(vals):
            raise ValueError("'dates' and 'vals' must be the same length.")
        if not all(isinstance(d, str) for d in dates):
            raise ValueError("All dates must be strings in 'yyyy-mm-dd' format.")

        # Grouping
        year_labels = [d[:4] for d in dates]
        ops = _ListOps()
        year_groups = ops.group_consecutive_duplicates_ordered(year_labels)
        return ops.group_by_sublist_len(vals, year_groups)

    def group_year_single(self,dates,vals,year):
        """Group the values of an input list ('vals') for a single specified year ('year') using a corresponding list of dates ('dates') and return a new list.

        Parameters:

        dates (list): the list of dates used to aggregate the input list of values. List of dates must be in a 'yyyy-mm-dd' format. Any character may be used to separate date parts. List of dates must be the same length as the list of values. It is expected that dates are ordered sequentially.

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.

        year (integer): the single year of data to group. Must be in a 'yyyy' format."""
        if not isinstance(dates, list) or not isinstance(vals, list):
            raise ValueError("Both 'dates' and 'vals' must be lists.")
        if len(dates) != len(vals):
            raise ValueError("'dates' and 'vals' must be the same length.")
        if not all(isinstance(d, str) for d in dates):
            raise ValueError("All dates must be strings in 'yyyy-mm-dd' format.")
        if not isinstance(year, int):
            raise ValueError("'year' must be an integer in 'yyyy' format.")

        # Grouping
        output = [v for d, v in zip(dates, vals) if int(d[:4]) == year]
        return output

    def group_by_month(self,dates,vals):
        """Group the values from an input list ('vals') using a corresponding list of dates ('dates') into lists of values within each month in the date range and return a nested list of monthly sublists.

        Parameters:

        dates (list): the list of dates used to aggregate the input list of values. List of dates must be in a 'yyyy-mm-dd' format. Any character may be used to separate date parts. List of dates must be the same length as the list of values. It is expected that dates are ordered sequentially.

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(dates, list) or not isinstance(vals, list):
            raise ValueError("Both 'dates' and 'vals' must be lists.")
        if len(dates) != len(vals):
            raise ValueError("'dates' and 'vals' must be the same length.")
        if not all(isinstance(d, str) for d in dates):
            raise ValueError("All dates must be strings in 'yyyy-mm-dd' format.")

        # Grouping
        month_labels = [d[:7] for d in dates]  # yyyy-mm
        ops = _ListOps()
        month_groups = ops.group_consecutive_duplicates_ordered(month_labels)
        return ops.group_by_sublist_len(vals, month_groups)
  
    def group_by_season(self,dates,vals,months):
        """Group the values from an input list ('vals') using a corresponding list of dates ('dates') into lists of values within the specified season ('months') in the date range and return a nested list of season sublists.

        Parameters:

        dates (list): the list of dates used to aggregate the input list of values. List of dates must be in a 'yyyy-mm-dd' format. Any character may be used to separate date parts. List of dates must be the same length as the list of values. It is expected that dates are ordered sequentially.

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation.
        
        months (list): a list of integers corresponding to the months within a particular season. For example, if grouping by the winter season, [1,2,12] would correspond to the months of January, February, and December, respectively. Month values should be given in ascending order."""
        if not isinstance(dates, list) or not isinstance(vals, list):
            raise ValueError("Both 'dates' and 'vals' must be lists.")
        if len(dates) != len(vals):
            raise ValueError("'dates' and 'vals' must be the same length.")
        if not all(isinstance(d, str) for d in dates):
            raise ValueError("All dates must be strings in 'yyyy-mm-dd' format.")
        if not isinstance(months, list) or not all(isinstance(m, int) for m in months):
            raise ValueError("'months' must be a list of integers.")

        # Filter by season months (single pass)
        season_dates, season_vals = [], []
        for d, v in zip(dates, vals):
            if int(d[5:7]) in months:
                season_dates.append(d)
                season_vals.append(v)

        # Group by year within the season
        return self.group_by_year(season_dates, season_vals)
   
    def group_by_day(self,dates,vals):
        """Group the values from an input list ('vals') using a corresponding list of dates ('dates') into lists of values within each day in the date range and return a nested list of daily sublists.

        Parameters:

        dates (list): the list of dates used to aggregate the input list of values. List of dates must be in a 'yyyy-mm-dd' format. Any character may be used to separate date parts. List of dates must be the same length as the list of values. It is expected that dates are ordered sequentially.

        vals (list): the list of values the operation will be performed on. The input list will not be altered by the operation."""
        if not isinstance(dates, list) or not isinstance(vals, list):
            raise ValueError("Both 'dates' and 'vals' must be lists.")
        if len(dates) != len(vals):
            raise ValueError("'dates' and 'vals' must be the same length.")
        if not all(isinstance(d, str) for d in dates):
            raise ValueError("All dates must be strings in 'yyyy-mm-dd' format.")

        # Grouping
        ops = _ListOps()
        day_groups = ops.group_consecutive_duplicates_ordered(dates)
        return ops.group_by_sublist_len(vals, day_groups)
  
class _Plot():
    """A set of functions for visualizing lists of data with common plot types."""
    # ============================================================
    # PRIVATE METHODS AND CLASS ATTRIBUTES AND FUNCTION DISPLAY
    # ============================================================
    def __init__(self):
        self._total_funcs = len([f for f in dir(__class__) if not f.startswith('__')])

    def __repr__(self):
        """If the toolbox is printed, display a message."""
        return ".plot: a set of functions for visualizing lists of data with common plot types."
    
    def _normalize_html_filename(self,filename):
        """Ensure the output filename of plots ends with .html.
        
        Parameters:
        
        filename (string): the filename to normalize."""
        directory = os.path.dirname(filename)
        base = os.path.basename(filename)
        root, ext = os.path.splitext(base)

        # Normalize to .html
        if ext.lower() != ".html":
            base = root + ".html"

        return os.path.join(directory, base)

    def _compute_categorical_ticks(self,n, num_ticks):
        """Compute the positions of categorical tick labels to sample.
        
        Parameters:
        
        n (integer): the length of the list to sample.
        
        num_ticks (integer): the number of ticks to derive samples for."""
        # Always include first and last
        num_ticks = num_ticks + 1
        positions = [0]

        if num_ticks > 2:
            k = num_ticks - 2
            step = (n - 2) / (k + 1)

            for i in range(1, k + 1):
                pos = int(round(1 + i * step))
                positions.append(pos)

        positions.append(n - 1)
        return positions

    def display_toolbox_functions(self):
        """Display a list of all available functions within this toolbox."""
        print(f"Number of {__class__.__name__[1:]} functions: {len([f for f in dir(__class__) if not f.startswith('__')])}")
        for f in [f for f in dir(__class__) if not f.startswith("__")]:
            print(f)

    # ============================================================
    # CONTINUOUS PLOTS
    # ============================================================

    def scatter_plot(self,x_var, y_vars,title="Scatter Plot",x_label="Independent Variable",y_label="Dependent Variable",num_ticks=5,point_size=3,point_colour=None, point_style=None,trendline=False,trend_colour=None,legend=False,point_labels=None,tooltip_threshold=100,margin=60, filename="scatter_plot.html",open_browser=False):
        """Create a scatter plot for a single independent variable ('x_var') and one or more dependent variables ('y_vars') and configure plot elements.

        Parameters:

        x_var (list): the list of values representing the independent variable to be visualized. Values may be numerical (integer or floating point) or categorical strings such as dates.

        y_vars (list): the list of values representing the dependent variable to be visualized, or a nested list if there are multiple dependent variables to be visualized. Values must be numerical (integer or floating point). Each series must match the length of x_var.

        title (string): the title of the scatter plot. 

        x_label (string): the label of the x axis.

        y_label (string): the label of the y axis. 

        num_ticks (integer): the number of tick marks to draw on each axis.

        point_size (integer): the radius of points on the plot (in pixels).

        point_colour (string, list): if a single dependent variable is given, this should be a single string representing the colour of the points. If multiple dependent variables are given, this should be a list of strings corresponding to the colour of each variable. If no colours are specified or the number of specified colours does not equal the number of dependent variables, random colours will be generated. Colours may be specified as a common colour name or hex code.

        point_style (string, list): the shape of points. May be specified as "circle", "diamond", "square", or "triangle". If a single dependent variable is given, this should be a single string representing the style of the points. If multiple dependent variables are given, this should be a list of strings corresponding to the style of each variable. If no styles are specified or the number of specified styles does not equal the number of dependent variables, random styles will be generated.

        trendline (Boolean): flag to indicate if trend lines for each dependent variable will be calculated and displayed.

        trend_colour (string, list): if a single dependent variable is given, this should be a single string representing the colour of the trend line. If multiple dependent variables are given, this should be a list of strings corresponding to the colour of each trend line. If no colours are specified or the number of specified colours does not equal the number of dependent variables, random colours will be generated. Colours may be specified as a common colour name or hex code.
        
        legend (Boolean): flag to indicate if a legend will be displayed.

        point_labels (string, list): if a single dependent variable is given, this should be a single string representing the label of the points. If multiple dependent variables are given, this should be a list of strings representing the labels of each variable. If no labels are given, generic labels will be generated. A legend must be created to view point labels.

        tooltip_threshold (integer): the maximum number of points per series for which explicit hover-based tooltips are included.
        
        margin (integer): the margin size (in pixels) around the plot area. 
        
        filename (string): the name of the output file. Include the output file path in the string or prepend tt.get_out_dir() to the filename.

        open_browser (Boolean): flag to indicate if the generated file should be automatically opened in the default web browser."""
        if not isinstance(x_var, list) or not x_var:
            raise ValueError("'x_var' must be a non-empty list of values.")
        x_is_cat = False
        if all(isinstance(v,str) for v in x_var):
            cat_ticks = [x_var[i] for i in self._compute_categorical_ticks(len(x_var),num_ticks)]
            x_var = [i for i in range(len(x_var))]
            x_is_cat = True
        else:
            if not all(isinstance(v, (int, float)) for v in x_var):
                raise ValueError("All values in 'x_var' must be numeric (int or float).")
        if not isinstance(y_vars, list):
            raise ValueError("'y_vars' must be a list of numeric values or a list of lists of numeric values.")
        is_flat = all(isinstance(v, (int, float)) for v in y_vars)
        is_nested = all(isinstance(v, list) for v in y_vars)
        if not is_flat and not is_nested:
            raise ValueError("'y_vars' must be a list of numeric values or a list of lists of numeric values.")
        if is_flat:
            if len(y_vars) != len(x_var):
                raise ValueError("The dependent variable series in 'y_vars' must match the length of 'x_var'.")
            num_series = 1
        else:
            if any(len(series) != len(x_var) for series in y_vars):
                raise ValueError("Each dependent variable series in 'y_vars' must match the length of 'x_var'.")
            for series in y_vars:
                if not all(isinstance(v, (int, float)) for v in series):
                    raise ValueError("All values in each dependent variable series must be numeric (int or float).")
            num_series = len(y_vars)
        if not isinstance(title, str):
            raise ValueError("'title' must be a string.")
        if not isinstance(x_label, str):
            raise ValueError("'x_label' must be a string.")
        if not isinstance(y_label, str):
            raise ValueError("'y_label' must be a string.")
        if not isinstance(num_ticks, int) or num_ticks <= 0:
            raise ValueError("'num_ticks' must be a positive integer.")
        if not isinstance(point_size, int) or point_size <= 0:
            raise ValueError("'point_size' must be a positive integer.")
        if point_colour is not None:
            if isinstance(point_colour, list):
                if len(point_colour) != num_series:
                    raise ValueError("If 'point_colour' is a list, it must match the number of dependent variable series.")
                if not all(isinstance(c, str) for c in point_colour):
                    raise ValueError("All entries in 'point_colour' must be strings representing colours.")
            elif not isinstance(point_colour, str):
                raise ValueError("'point_colour' must be a string or a list of strings.")
        valid_styles = {"circle", "square", "diamond", "triangle"}
        if point_style is not None:
            if isinstance(point_style, list):
                if len(point_style) != num_series:
                    raise ValueError("If 'point_style' is a list, it must match the number of dependent variable series.")
                if not all(isinstance(s, str) and s.lower() in valid_styles for s in point_style):
                    raise ValueError("Each point style must be one of: 'circle', 'square', 'diamond', 'triangle'.")
            elif isinstance(point_style, str):
                if point_style.lower() not in valid_styles:
                    raise ValueError("'point_style' must be one of: 'circle', 'square', 'diamond', 'triangle'.")
            else:
                raise ValueError("'point_style' must be a string or a list of strings.")
        if not isinstance(trendline, bool):
            raise ValueError("'trendline' must be a Boolean value.")
        if trend_colour is not None:
            if isinstance(trend_colour, list):
                if len(trend_colour) != num_series:
                    raise ValueError("If 'trend_colour' is a list, it must match the number of dependent variable series.")
                if not all(isinstance(c, str) for c in trend_colour):
                    raise ValueError("All entries in 'trend_colour' must be strings representing colours.")
            elif not isinstance(trend_colour, str):
                raise ValueError("'trend_colour' must be a string or a list of strings.")
        if not isinstance(legend, bool):
            raise ValueError("'legend' must be a Boolean value.")
        if point_labels is not None:
            if isinstance(point_labels, list):
                if len(point_labels) != num_series:
                    raise ValueError("If 'point_labels' is a list, it must match the number of dependent variable series.")
                if not all(isinstance(lbl, str) for lbl in point_labels):
                    raise ValueError("All entries in 'point_labels' must be strings.")
            elif not isinstance(point_labels, str):
                raise ValueError("'point_labels' must be a string or a list of strings.")
        if not isinstance(tooltip_threshold, int) or tooltip_threshold < 0:
            raise ValueError("'tooltip_threshold' must be a non-negative integer.")
        if not isinstance(margin, int) or margin < 0:
            raise ValueError("'margin' must be a non-negative integer.")
        if not isinstance(filename, str):
            raise ValueError("'filename' must be a string.")
        filename = self._normalize_html_filename(filename)
        if not isinstance(open_browser, bool):
            raise ValueError("'open_browser' must be a Boolean value.")
        
        # Normalize y_vars to a list of lists
        if isinstance(y_vars[0], (list, tuple)):
            series_list = y_vars
        else:
            series_list = [y_vars]

        num_series = len(series_list)

        # Strict enforcement: all series must match x_var length
        for idx, series in enumerate(series_list):
            if len(series) != len(x_var):
                raise ValueError(
                    f"Length mismatch: series {idx} has {len(series)} points, "
                    f"but x_var has {len(x_var)} points."
                )

        # Dimensions and ranges
        base_width, base_height = 800, 600

        xmin, xmax = min(x_var), max(x_var)
        ymin = min(min(series) for series in series_list)
        ymax = max(max(series) for series in series_list)

        if xmax == xmin: xmax = xmin + 1
        if ymax == ymin: ymax = ymin + 1

        # Estimate tick label width for dynamic Y-label offset
        tick_labels = [f"{ymin + i * (ymax - ymin) / num_ticks:.1f}" for i in range(num_ticks + 1)]
        max_label_len = max(len(lbl) for lbl in tick_labels)
        font_size = 12  # matches tick label font size
        est_label_width = max_label_len * font_size * 0.6
        y_label_offset = margin - est_label_width - 15  # 15 px padding

        # Reserve extra canvas space so label + legend are inside SVG
        extra_left = max(0, margin - y_label_offset)  # ensures leftward label fits
        extra_right = 180 if legend else 0            # ensures legend fits

        width = base_width + margin*2 + extra_left + extra_right
        height = base_height + margin*2

        # Plot area
        plot_left   = margin + extra_left
        plot_right  = width - margin - extra_right
        plot_top    = margin
        plot_bottom = height - margin

        # Scales
        xscale = (plot_right - plot_left) / (xmax - xmin)
        yscale = (plot_bottom - plot_top) / (ymax - ymin)

        sx = lambda x: plot_left + (x - xmin) * xscale
        sy = lambda y: plot_bottom - (y - ymin) * yscale

        # Style normalization
        def normalize_style(user_input, generator, label):
            if user_input is None:
                return [generator() for _ in range(num_series)]
            elif isinstance(user_input, str):
                return [user_input] * num_series
            elif isinstance(user_input, (list, tuple)):
                return [user_input[i % len(user_input)] for i in range(num_series)]
            else:
                raise ValueError(f"Invalid {label} input")

        def random_colour_palette(n_series):
            base_step = 360 / n_series
            offset = random.randint(0, 359)  # random starting hue
            return [
                f"hsl({int((offset + i * base_step) % 360)},70%,50%)"
                for i in range(n_series)
            ]

        def random_style():
            return random.choice(["circle", "square", "diamond", "triangle"])

        # Normalize colours and styles
        if point_colour is None:
            point_colours = random_colour_palette(num_series)
        else:
            point_colours = normalize_style(point_colour, None, "point_colour")

        point_styles = normalize_style(point_style, random_style, "point_style")

        if trend_colour is None:
            trend_colours = point_colours
        else:
            trend_colours = normalize_style(trend_colour, None, "trend_colour")


        # Normalize series labels
        if point_labels is None:
            point_labels = [f"Series {i+1}" for i in range(num_series)]
        elif isinstance(point_labels,str):
            point_labels = [point_labels]
        else:
            point_labels = [point_labels[i] if i < len(point_labels) else f"Series {i+1}"
                            for i in range(num_series)]

        # Begin SVG (fixed)
        buf = io.StringIO()
        buf.write(
            f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}" '
            f'preserveAspectRatio="xMidYMid meet">\n'
        )
        buf.write('  <rect width="100%" height="100%" fill="white"/>\n')

        buf.write(f'  <title>{title}</title>\n')
        buf.write(f'  <desc>Scatterplot with {num_series} series.</desc>\n')

        # Definitions for point styles
        buf.write('  <defs>\n')
        for idx, (style, colour) in enumerate(zip(point_styles, point_colours)):
            if style == "circle":
                shape_def = f'<circle id="dot{idx}" r="{point_size}" fill="{colour}"/>'
            elif style == "square":
                half = point_size
                shape_def = f'<rect id="dot{idx}" x="-{half}" y="-{half}" width="{2*half}" height="{2*half}" fill="{colour}"/>'
            elif style == "diamond":
                s = point_size
                shape_def = f'<polygon id="dot{idx}" points="0,-{s} {s},0 0,{s} -{s},0" fill="{colour}"/>'
            elif style == "triangle":
                s = point_size
                shape_def = f'<polygon id="dot{idx}" points="0,-{s} {s},{s} -{s},{s}" fill="{colour}"/>'
            buf.write(f'    {shape_def}\n')
        buf.write('  </defs>\n')

        # Axes using plot bounds (not raw margin)
        buf.write(f'  <line x1="{plot_left}" y1="{plot_bottom}" x2="{plot_right}" y2="{plot_bottom}" stroke="black"/>\n')  # X axis
        buf.write(f'  <line x1="{plot_left}" y1="{plot_top}"    x2="{plot_left}"  y2="{plot_bottom}" stroke="black"/>\n')  # Y axis

        # Title and X label
        buf.write(f'  <text x="{(plot_left + plot_right)/2}" y="{height-15}" text-anchor="middle" font-size="16">{x_label}</text>\n')
        buf.write(f'  <text x="{width/2}" y="30" text-anchor="middle" font-size="20" font-weight="bold">{title}</text>\n')

        # Ticks and tick labels using plot area
        for i in range(num_ticks + 1):
            # X ticks
            tx = plot_left + i * (plot_right - plot_left) / num_ticks
            if x_is_cat:
                valx = cat_ticks[i]
                buf.write(f'  <text x="{tx}" y="{plot_bottom+20}" text-anchor="middle" font-size="12">{valx}</text>\n')
            else:
                valx = xmin + i * (xmax - xmin) / num_ticks
                buf.write(f'  <text x="{tx}" y="{plot_bottom+20}" text-anchor="middle" font-size="12">{valx:.1f}</text>\n')
            # valx = xmin + i * (xmax - xmin) / num_ticks
            buf.write(f'  <line x1="{tx}" y1="{plot_bottom}" x2="{tx}" y2="{plot_bottom+5}" stroke="black"/>\n')

            # Y ticks
            ty = plot_bottom - i * (plot_bottom - plot_top) / num_ticks
            valy = ymin + i * (ymax - ymin) / num_ticks
            buf.write(f'  <line x1="{plot_left}" y1="{ty}" x2="{plot_left-5}" y2="{ty}" stroke="black"/>\n')
            buf.write(f'  <text x="{plot_left-10}" y="{ty+4}" text-anchor="end" font-size="12">{valy:.1f}</text>\n')

        # Rotated Y-axis label: keep inside extra_left padding
        y_label_x = plot_left - (extra_left - 5)  # 5px inset from left padding
        y_label_y = (plot_top + plot_bottom) / 2
        buf.write(
            f'  <text x="{y_label_x}" y="{y_label_y}" text-anchor="middle" font-size="16" '
            f'transform="rotate(-90,{y_label_x},{y_label_y})">{y_label}</text>\n'
        )

        # Legend
        if legend:
            legend_font = 12
            swatch_size = 14
            swatch_gap = 10
            legend_row_gap = 6

            # Estimate text widths for labels (include "trend" if enabled)
            est_text_widths = []
            for lbl in point_labels:
                base_width = len(str(lbl)) * 7
                if trendline:
                    trend_width = len(str(lbl) + " trend") * 7
                    est_text_widths.append(max(base_width, trend_width))
                else:
                    est_text_widths.append(base_width)

            legend_inner_width = max(swatch_size + swatch_gap + tw for tw in est_text_widths)

            # Ensure the "Legend" title always fits
            min_title_width = len("Legend") * 8 + 20
            legend_width = max(legend_inner_width + 20, min_title_width)

            # Height accounts for series entries (and trendline rows if enabled)
            legend_height = len(point_labels) * (swatch_size + legend_row_gap) * (2 if trendline else 1) + 30

            legend_x = plot_right + 40
            legend_y = plot_top

            # Legend frame with rounded corners
            buf.write(
                f"<rect x='{legend_x}' y='{legend_y}' width='{legend_width}' height='{legend_height}' "
                f"fill='white' stroke='black' rx='6' ry='6'/>\n"
            )
            buf.write(
                f"<text x='{legend_x + 10}' y='{legend_y + 18}' font-size='14' text-anchor='start'>Legend</text>\n"
            )

            # Entries
            for idx, (style, colour, label) in enumerate(zip(point_styles, point_colours, point_labels)):
                entry_y = legend_y + 30 + legend_row_gap + idx * (swatch_size + legend_row_gap) * (2 if trendline else 1)

                # Draw symbol for series
                if style == "circle":
                    buf.write(f"<circle cx='{legend_x + 10}' cy='{entry_y}' r='{point_size}' fill='{colour}' stroke='black'/>\n")
                elif style == "square":
                    half = point_size
                    buf.write(f"<rect x='{legend_x + 10 - half}' y='{entry_y - half}' width='{2*half}' height='{2*half}' fill='{colour}' stroke='black'/>\n")
                elif style == "diamond":
                    s = point_size
                    buf.write(f"<polygon points='{legend_x+10},{entry_y-s} {legend_x+10+s},{entry_y} {legend_x+10},{entry_y+s} {legend_x+10-s},{entry_y}' fill='{colour}' stroke='black'/>\n")
                elif style == "triangle":
                    s = point_size
                    buf.write(f"<polygon points='{legend_x+10},{entry_y-s} {legend_x+10+s},{entry_y+s} {legend_x+10-s},{entry_y+s}' fill='{colour}' stroke='black'/>\n")

                buf.write(f"<text x='{legend_x + 10 + swatch_size + swatch_gap}' y='{entry_y+4}' font-size='{legend_font}' text-anchor='start'>{label}</text>\n")

                # Trend line entry
                if trendline:
                    y_trend = entry_y + swatch_size + legend_row_gap
                    buf.write(f"<line x1='{legend_x + 5}' y1='{y_trend}' x2='{legend_x + 25}' y2='{y_trend}' stroke='{trend_colours[idx]}' stroke-width='2'/>\n")
                    buf.write(f"<text x='{legend_x + 10 + swatch_size + swatch_gap}' y='{y_trend+4}' font-size='{legend_font}' text-anchor='start'>{label} trend</text>\n")

            buf.write('  </g>\n')

        # Series points (use sx/sy consistently)
        for idx, series in enumerate(series_list):
            buf.write(f'  <g id="series{idx}">\n')

            include_tooltips = len(series) <= tooltip_threshold
            style = point_styles[idx]
            colour = point_colours[idx]

            if include_tooltips:
                for x, y in zip(x_var, series):
                    cx, cy = sx(x), sy(y)
                    if style == "circle":
                        buf.write(f'    <circle cx="{cx}" cy="{cy}" r="{point_size}" fill="{colour}"><title>({x:.2f},{y:.2f})</title></circle>\n')
                    elif style == "square":
                        half = point_size
                        buf.write(f'    <rect x="{cx-half}" y="{cy-half}" width="{2*half}" height="{2*half}" fill="{colour}"><title>({x:.2f},{y:.2f})</title></rect>\n')
                    elif style == "diamond":
                        s = point_size
                        buf.write(f'    <polygon points="{cx},{cy-s} {cx+s},{cy} {cx},{cy+s} {cx-s},{cy}" fill="{colour}"><title>({x:.2f},{y:.2f})</title></polygon>\n')
                    elif style == "triangle":
                        s = point_size
                        buf.write(f'    <polygon points="{cx},{cy-s} {cx+s},{cy+s} {cx-s},{cy+s}" fill="{colour}"><title>({x:.2f},{y:.2f})</title></polygon>\n')
            else:
                # Compact <use> references (no tooltips)
                point_lines = "".join(
                    f'    <use href="#dot{idx}" x="{sx(x)}" y="{sy(y)}"/>\n'
                    for x, y in zip(x_var, series)
                )
                buf.write(point_lines)

            buf.write('  </g>\n')

            # Trend line with tooltip (use sx/sy and plot bounds)
            if trendline:
                n = len(x_var)
                sum_x = sum(x_var)
                sum_y = sum(series)
                sum_xy = sum(x * y for x, y in zip(x_var, series))
                sum_x2 = sum(x * x for x in x_var)

                denom = (n * sum_x2 - sum_x ** 2)
                slope = (n * sum_xy - sum_x * sum_y) / denom if denom != 0 else 0.0
                intercept = (sum_y - slope * sum_x) / n

                x1_data, y1_data = xmin, slope * xmin + intercept
                x2_data, y2_data = xmax, slope * xmax + intercept

                def clip_to_y(xd, yd):
                    if yd < ymin:
                        xd = (ymin - intercept) / slope if slope != 0 else xd
                        yd = ymin
                    elif yd > ymax:
                        xd = (ymax - intercept) / slope if slope != 0 else xd
                        yd = ymax
                    return xd, yd

                x1_data, y1_data = clip_to_y(x1_data, y1_data)
                x2_data, y2_data = clip_to_y(x2_data, y2_data)

                x1_data = max(min(x1_data, xmax), xmin)
                x2_data = max(min(x2_data, xmax), xmin)

                x1 = sx(x1_data)
                y1 = sy(y1_data)
                x2 = sx(x2_data)
                y2 = sy(y2_data)

                buf.write(
                    f'  <line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" '
                    f'stroke="{trend_colours[idx]}" stroke-width="2">'
                    f'<title>y = {slope:.3f}x + {intercept:.3f}</title></line>\n'
                )
        # End SVG
        buf.write('</svg>\n')

        # Also write HTML wrapper with save button
        html_filename = filename
        html_content = f"""<!DOCTYPE html>
        <html>
        <head>
        <meta charset="utf-8">
        <title>{title}</title>
        <style>
            html, body {{
            height: 100%;
            margin: 0;
            }}
            .wrap {{
            width: 100vw;
            height: 100vh;
            display: flex;
            justify-content: center;
            align-items: center;
            }}
            svg {{
            max-width: 100vw;
            max-height: 100vh;
            width: 100%;
            height: auto;
            display: block;
            background: white;
            font-family: Arial, Helvetica, sans-serif;
            }}
            #saveBtn {{
            position: fixed;
            top: 10px;
            right: 10px;
            padding: 6px 12px;
            font-size: 14px;
            cursor: pointer;
            }}
        </style>
        </head>
        <body>
        <div class="wrap">
            {buf.getvalue()}
        </div>
        <button id="saveBtn" onclick="downloadPNG()">Save as PNG</button>
        <script>
            function downloadPNG() {{
            const svg = document.querySelector("svg");
            const serializer = new XMLSerializer();
            const svgStr = serializer.serializeToString(svg);

            // Use full viewBox dimensions if present
            let width, height;
            if (svg.viewBox && svg.viewBox.baseVal) {{
                width = svg.viewBox.baseVal.width;
                height = svg.viewBox.baseVal.height;
            }} else {{
                width = svg.width.baseVal.value;
                height = svg.height.baseVal.value;
            }}

            const canvas = document.createElement("canvas");
            canvas.width = width;
            canvas.height = height;
            const ctx = canvas.getContext("2d");

            const img = new Image();
            img.onload = function() {{
                ctx.drawImage(img, 0, 0, width, height);
                const link = document.createElement("a");
                // Use plot title instead of file path
                const safeTitle = "{title}".replace(/[^a-z0-9]+/gi, "_").toLowerCase();
                link.download = safeTitle + ".png";
                link.href = canvas.toDataURL("image/png");
                link.click();
            }};
            img.src = "data:image/svg+xml;base64," + btoa(unescape(encodeURIComponent(svgStr)));
            }}
        </script>
        </body>
        </html>"""

        # Write file
        with open(html_filename, "w") as f:
            f.write(html_content)

        # Open in browser
        if open_browser:
            import webbrowser
            webbrowser.open(html_filename)

    def line_plot(self,x_var, y_vars,title="Line Plot",x_label="Independent Variable", y_label="Dependent Variable",num_ticks=5,line_weight=2,line_colours=None,trendline=False,trend_colours=None,legend=False,line_labels=None,margin=60,filename="line_plot.html",open_browser=False):
        """Create a line plot for a single independent variable ('x_var') and one or more dependent variables ('y_vars')  and configure plot elements.

        Parameters:

        x_var (list): the list of values representing the independent variable to be visualized. Values may be numerical (integer or floating point) or categorical strings such as dates.

        y_vars (list): the list of values representing the dependent variable to be visualized, or a nested list if there are multiple dependent variables to be visualized. Values must be numerical (integer or floating point). Each series must match the length of x_var.

        title (string): the title of the line plot.  

        x_label (string): the label of the x axis.

        y_label (string): the label of the y axis.  

        num_ticks (integer): the number of tick marks to draw on each axis.

        line_weight (integer): the stroke width of the lines on the plot (in pixels).

        line_colours (string, list): if a single dependent variable is given, this should be a single string representing the colour of the line. If multiple dependent variables are given, this should be a list of strings corresponding to the colour of each variable. If no colours are specified or the number of specified colours does not equal the number of dependent variables, random colours will be generated. Colours may be specified as a common colour name or hex code.

        trendline (Boolean): flag to indicate if trend lines for each dependent variable will be calculated and displayed.

        trend_colours (string, list): if a single dependent variable is given, this should be a single string representing the colour of the trend line. If multiple dependent variables are given, this should be a list of strings corresponding to the colour of each trend line. If no colours are specified or the number of specified colours does not equal the number of dependent variables, random colours will be generated. Colours may be specified as a common colour name or hex code.
        
        legend (Boolean): flag to indicate if a legend will be displayed.

        line_labels (string, list): the labels of each dependent variable series. If no labels are given, generic labels ("Series 1", "Series 2", ...) will be generated. A legend must be created to view series labels.
        
        margin (integer): the margin size (in pixels) around the plot area.  
        
        filename (string): the name of the output file. Include the output file path in the string or prepend tt.get_out_dir() to the filename.

        open_browser (Boolean): flag to indicate if the generated file should be automatically opened in the default web browser."""
        if not isinstance(x_var, list) or not x_var:
            raise ValueError("'x_var' must be a non-empty list of values.")
        x_is_cat = False
        if all(isinstance(v,str) for v in x_var):
            cat_ticks = [x_var[i] for i in self._compute_categorical_ticks(len(x_var),num_ticks)]
            x_var = [i for i in range(len(x_var))]
            x_is_cat = True
        else:
            if not all(isinstance(v, (int, float)) for v in x_var):
                raise ValueError("All values in 'x_var' must be numeric (int or float).")
        if not isinstance(y_vars, list):
            raise ValueError("'y_vars' must be a list of numeric values or a list of lists of numeric values.")
        is_flat = all(isinstance(v, (int, float)) for v in y_vars)
        is_nested = all(isinstance(v, list) for v in y_vars)
        if not is_flat and not is_nested:
            raise ValueError("'y_vars' must be a list of numeric values or a list of lists of numeric values.")
        if is_flat:
            if len(y_vars) != len(x_var):
                raise ValueError("The dependent variable series in 'y_vars' must match the length of 'x_var'.")
            num_series = 1
        else:
            if any(len(series) != len(x_var) for series in y_vars):
                raise ValueError("Each dependent variable series in 'y_vars' must match the length of 'x_var'.")
            for series in y_vars:
                if not all(isinstance(v, (int, float)) for v in series):
                    raise ValueError("All values in each dependent variable series must be numeric (int or float).")
            num_series = len(y_vars)
        if not isinstance(title, str):
            raise ValueError("'title' must be a string.")
        if not isinstance(x_label, str):
            raise ValueError("'x_label' must be a string.")
        if not isinstance(y_label, str):
            raise ValueError("'y_label' must be a string.")
        if not isinstance(num_ticks, int) or num_ticks <= 0:
            raise ValueError("'num_ticks' must be a positive integer.")
        if not isinstance(line_weight, int) or line_weight <= 0:
            raise ValueError("'line_weight' must be a positive integer.")
        if line_colours is not None:
            if isinstance(line_colours, list):
                if len(line_colours) != num_series:
                    raise ValueError("If 'line_colours' is a list, it must match the number of dependent variable series.")
                if not all(isinstance(c, str) for c in line_colours):
                    raise ValueError("All entries in 'line_colours' must be strings representing colours.")
            elif not isinstance(line_colours, str):
                raise ValueError("'line_colours' must be a string or a list of strings.")
        if not isinstance(trendline, bool):
            raise ValueError("'trendline' must be a Boolean value.")
        if trend_colours is not None:
            if isinstance(trend_colours, list):
                if len(trend_colours) != num_series:
                    raise ValueError("If 'trend_colours' is a list, it must match the number of dependent variable series.")
                if not all(isinstance(c, str) for c in trend_colours):
                    raise ValueError("All entries in 'trend_colours' must be strings representing colours.")
            elif not isinstance(trend_colours, str):
                raise ValueError("'trend_colours' must be a string or a list of strings.")
        if not isinstance(legend, bool):
            raise ValueError("'legend' must be a Boolean value.")
        if line_labels is not None:
            if not isinstance(line_labels, (str, list)):
                raise ValueError("'line_labels' must be a string or a list of strings.")
            if len(line_labels) != num_series and not isinstance(line_labels,str):
                raise ValueError("Length of 'line_labels' must match the number of dependent variable series.")
            if not all(isinstance(lbl, str) for lbl in line_labels):
                raise ValueError("All entries in 'line_labels' must be strings.")
        if not isinstance(margin, int) or margin < 0:
            raise ValueError("'margin' must be a non-negative integer.")
        if not isinstance(filename, str):
            raise ValueError("'filename' must be a string.")
        filename = self._normalize_html_filename(filename)
        if not isinstance(open_browser, bool):
            raise ValueError("'open_browser' must be a Boolean value.")
        
        # sort values to ensure correct drawing
        order = sorted(range(len(x_var)), key=lambda i: x_var[i])
        x_var = [x_var[i] for i in order]
        
        if is_flat:
            y_vars = [y_vars[i] for i in order]
        else:
            y_vars = [
                [series[i] for i in order]
                for series in y_vars
            ]

        # Normalize y_vars to list of series
        if not isinstance(y_vars[0], (list, tuple)):
            series_list = [y_vars]
        else:
            series_list = y_vars
        num_series = len(series_list)

        # Normalize colours
        if line_colours is None:
            base_step = 360 / num_series
            offset = random.randint(0, 359)  # random starting hue
            line_colours = [
                f"hsl({int((offset + i * base_step) % 360)},70%,40%)"
                for i in range(num_series)
            ]

        # Match trend colours to line colours unless user supplied
        if trend_colours is None:
            trend_colours = line_colours

        if line_labels is None:
            line_labels = [f"Series {i+1}" for i in range(num_series)]
        elif isinstance(line_labels,str):
            line_labels = [line_labels]

        # Data ranges
        xmin, xmax = min(x_var), max(x_var)
        ymin = min(min(series) for series in series_list)
        ymax = max(max(series) for series in series_list)

        # Estimate tick label width for dynamic Y-label offset
        tick_labels = [f"{round(ymin + i * (ymax - ymin) / num_ticks, 2)}" for i in range(num_ticks + 1)]
        max_label_len = max(len(lbl) for lbl in tick_labels)
        font_size = 12  # matches tick label font size
        est_label_width = max_label_len * font_size * 0.6
        y_label_offset = margin - est_label_width - 15  # 15 px padding

        # Reserve extra canvas space so label fits
        extra_left = max(0, margin - y_label_offset)
        extra_left = est_label_width + 50
        extra_right = margin + (100 if legend else 0)
        extra_top = margin + 40
        extra_bottom = margin + 40

        width = 800 + extra_left + extra_right
        height = 600 + extra_top + extra_bottom

        plot_left = margin + extra_left
        plot_right = width - extra_right
        plot_top = extra_top
        plot_bottom = height - extra_bottom


        # Coordinate transforms
        def sx(x): return plot_left + (x - xmin) / (xmax - xmin) * (plot_right - plot_left)
        def sy(y): return plot_bottom - (y - ymin) / (ymax - ymin) * (plot_bottom - plot_top)

        # Begin SVG output
        buffer = io.StringIO()
        buffer.write(f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' "
                    f"viewBox='0 0 {width} {height}' style='font-family:sans-serif'>\n")
        
        buffer.write('  <rect width="100%" height="100%" fill="white"/>\n')

        # Title and description for accessibility
        buffer.write(f"<title>{title}</title>\n")
        buffer.write(f"<desc>Responsive line plot with {num_series} series</desc>\n")

        # Axes
        buffer.write(f"<line x1='{plot_left}' y1='{plot_top}' x2='{plot_left}' y2='{plot_bottom}' "
                    f"stroke='black'/>\n")
        buffer.write(f"<line x1='{plot_left}' y1='{plot_bottom}' x2='{plot_right}' y2='{plot_bottom}' "
                    f"stroke='black'/>\n")

        # Axis labels
        y_label_x = (plot_left - 10) - est_label_width  - 20 # 10px inset from left padding
        y_label_y = (plot_top + plot_bottom) / 2

        buffer.write(f"<text x='{(plot_left+plot_right)/2}' y='{height - margin}' "
                    f"text-anchor='middle'>{x_label}</text>\n")
        buffer.write(f"<text x='{y_label_x}' y='{y_label_y}' text-anchor='middle' "
                f"transform='rotate(-90 {y_label_x},{y_label_y})'>{y_label}</text>\n")

        # Title
        buffer.write(f"<text x='{(plot_left+plot_right)/2}' y='{margin}' "
                    f"text-anchor='middle' font-size='16'>{title}</text>\n")

        # Ticks and tick labels
        for i in range(num_ticks+1):
            # X ticks
            xt = xmin + i*(xmax-xmin)/num_ticks
            x_pos = sx(xt)
            if x_is_cat:
                xt = cat_ticks[i]
                buffer.write(f"<text x='{x_pos}' y='{plot_bottom+20}' text-anchor='middle'>{xt}</text>\n")
            else:
                buffer.write(f"<text x='{x_pos}' y='{plot_bottom+20}' text-anchor='middle'>{round(xt,2)}</text>\n")
            buffer.write(f"<line x1='{x_pos}' y1='{plot_bottom}' x2='{x_pos}' y2='{plot_bottom+5}' stroke='black'/>\n")

            # Y ticks
            yt = ymin + i*(ymax-ymin)/num_ticks
            y_pos = sy(yt)
            buffer.write(f"<line x1='{plot_left-5}' y1='{y_pos}' x2='{plot_left}' y2='{y_pos}' stroke='black'/>\n")
            buffer.write(f"<text x='{plot_left-10}' y='{y_pos+4}' text-anchor='end'>{round(yt,2)}</text>\n")

        #Legend
        if legend:
            legend_font = 12
            swatch_size = 20   # length of line segment
            swatch_gap = 10
            legend_row_gap = 6

            # Estimate text widths for labels (include "trend" if needed)
            est_text_widths = []
            for lbl in line_labels:
                base_width = len(str(lbl)) * 7
                trend_width = len(str(lbl) + " trend") * 7 if trendline else base_width
                est_text_widths.append(max(base_width, trend_width))

            legend_inner_width = max(swatch_size + swatch_gap + tw for tw in est_text_widths)

            # Ensure the "Legend" title always fits
            min_title_width = len("Legend") * 8 + 20
            legend_width = max(legend_inner_width + 20, min_title_width)

            # Height accounts for series entries (two rows per series if trend lines are shown)
            legend_height = len(line_labels) * (swatch_size + legend_row_gap) * (2 if trendline else 1) + 30

            legend_x = plot_right + 40
            legend_y = plot_top

            # Legend frame with rounded corners
            buffer.write(
                f"<rect x='{legend_x}' y='{legend_y}' width='{legend_width}' height='{legend_height}' "
                f"fill='white' stroke='black' rx='6' ry='6'/>\n"
            )
            buffer.write(
                f"<text x='{legend_x + 10}' y='{legend_y + 18}' font-size='14' text-anchor='start'>Legend</text>\n"
            )

            # Entries
            for i, label in enumerate(line_labels):
                entry_y = legend_y + 30 + legend_row_gap + i * (swatch_size + legend_row_gap) * (2 if trendline else 1)

                # Draw line swatch
                buffer.write(
                    f"<line x1='{legend_x + 10}' y1='{entry_y}' x2='{legend_x + 10 + swatch_size}' y2='{entry_y}' "
                    f"stroke='{line_colours[i]}' stroke-width='{line_weight}'/>\n"
                )
                buffer.write(
                    f"<text x='{legend_x + 10 + swatch_size + swatch_gap}' y='{entry_y+4}' "
                    f"font-size='{legend_font}' text-anchor='start'>{label}</text>\n"
                )

                # Trend line entry (second row)
                if trendline:
                    y_trend = entry_y + swatch_size + legend_row_gap
                    buffer.write(
                        f"<line x1='{legend_x + 10}' y1='{y_trend}' x2='{legend_x + 10 + swatch_size}' y2='{y_trend}' "
                        f"stroke='{trend_colours[i]}' stroke-dasharray='4' stroke-width='1'/>\n"
                    )
                    buffer.write(
                        f"<text x='{legend_x + 10 + swatch_size + swatch_gap}' y='{y_trend+4}' "
                        f"font-size='{legend_font}' text-anchor='start'>{label} trend</text>\n"
        )
        
        # Draw series lines and tooltips
        for i, series in enumerate(series_list):
            colour = line_colours[i]

            # Build path string for the line
            path_d = "M " + " ".join(f"{sx(x_var[j])},{sy(series[j])}" for j in range(len(x_var)))
            buffer.write(f"<path d='{path_d}' fill='none' stroke='{colour}' stroke-width='{line_weight}'/>\n")

            # Add minimal markers at start and end points with tooltips
            if len(series) > 0:
                # Start point
                x_start, y_start = sx(x_var[0]), sy(series[0])
                buffer.write(f"<circle cx='{x_start}' cy='{y_start}' r='1' fill='transparent' stroke='none'>"
                            f"<title>{line_labels[i]} start: ({x_var[0]}, {series[0]})</title></circle>\n")

                # End point
                x_end, y_end = sx(x_var[-1]), sy(series[-1])
                buffer.write(f"<circle cx='{x_end}' cy='{y_end}' r='1' fill='transparent' stroke='none'>"
                            f"<title>{line_labels[i]} end: ({x_var[-1]}, {series[-1]})</title></circle>\n")

            # Optional trendline
            if trendline:
                n = len(x_var)
                x_mean = sum(x_var)/n
                y_mean = sum(series)/n
                num = sum((x_var[j]-x_mean)*(series[j]-y_mean) for j in range(n))
                den = sum((x_var[j]-x_mean)**2 for j in range(n))
                slope = num/den if den != 0 else 0
                intercept = y_mean - slope*x_mean

                y_start = slope*xmin + intercept
                y_end = slope*xmax + intercept
                buffer.write(f"<line x1='{sx(xmin)}' y1='{sy(y_start)}' "
                            f"x2='{sx(xmax)}' y2='{sy(y_end)}' "
                            f"stroke='{trend_colours[i]}' stroke-dasharray='4' stroke-width='1'/>\n")

        # Close SVG
        buffer.write("</svg>\n")

        # Also write HTML wrapper with save button
        html_filename = filename
        html_content = f"""<!DOCTYPE html>
        <html>
        <head>
        <meta charset="utf-8">
        <title>{title}</title>
        <style>
            html, body {{
                height: 100%;
                margin: 0;
            }}
            .wrap {{
                width: 100vw;
                height: 100vh;
                display: flex;
                justify-content: center;
                align-items: center;
            }}
            svg {{
                max-width: 100vw;
                max-height: 100vh;
                width: 100%;
                height: auto;
                display: block;
                background: white;
                font-family: Arial, Helvetica, sans-serif;
            }}
            #saveBtn {{
                position: fixed;
                top: 10px;
                right: 10px;
                padding: 6px 12px;
                font-size: 14px;
                cursor: pointer;
            }}
        </style>
        </head>
        <body>
        <div class="wrap">
            {buffer.getvalue()}
        </div>
        <button id="saveBtn" onclick="downloadPNG()">Save as PNG</button>
        <script>
        function downloadPNG() {{
            const svg = document.querySelector("svg");
            const serializer = new XMLSerializer();
            const svgStr = serializer.serializeToString(svg);

            let width, height;
            if (svg.viewBox && svg.viewBox.baseVal) {{
                width = svg.viewBox.baseVal.width;
                height = svg.viewBox.baseVal.height;
            }} else {{
                width = svg.width.baseVal.value;
                height = svg.height.baseVal.value;
            }}

            const canvas = document.createElement("canvas");
            canvas.width = width;
            canvas.height = height;
            const ctx = canvas.getContext("2d");

            const img = new Image();
            img.onload = function() {{
                ctx.drawImage(img, 0, 0, width, height);
                const link = document.createElement("a");
                const safeTitle = "{title}".replace(/[^a-z0-9]+/gi, "_").toLowerCase();
                link.download = safeTitle + ".png";
                link.href = canvas.toDataURL("image/png");
                link.click();
            }};
            img.src = "data:image/svg+xml;base64," + btoa(unescape(encodeURIComponent(svgStr)));
        }}
        </script>
        </body>
        </html>"""

        with open(html_filename, "w", encoding="utf-8") as f:
            f.write(html_content)

        # Optionally open in browser (open HTML wrapper so button is visible)
        if open_browser:
            import webbrowser
            webbrowser.open(html_filename)

    # ============================================================
    # DISTRIBUTION PLOTS
    # ============================================================

    def histogram(self,variable,title="Histogram",x_label="Value",y_label="Frequency",num_bins=None,bar_colour=None,margin=60,filename="histogram.html",open_browser=False):
        """Create a histogram for a single data series ('variable') and configure plot elements.

        Parameters:

        variable (list): the list of numerical values to be visualized. Values must be integers or floating point.

        title (string): the title of the histogram.

        x_label (string): the label of the x axis.

        y_label (string): the label of the y axis.

        num_bins (integer, optional): the number of bins to divide the data into. If not specified, the optimal number of bins is calculated using Sturge's Formula.

        bar_colour (string, optional): the colour of the bars. If not specified, a random colour will be generated. Colours may be specified as a common colour name or hex code.

        margin (integer): the margin size (in pixels) around the plot area.

        filename (string): the name of the output file. Include the output file path in the string or prepend tt.get_out_dir() to the filename.

        open_browser (Boolean): flag to indicate if the generated file should be automatically opened in the default web browser."""
        if not isinstance(variable, list) or not variable:
            raise ValueError("'variable' must be a non-empty list of numeric values.")
        if not all(isinstance(v, (int, float)) for v in variable):
            raise ValueError("All values in 'variable' must be numeric (int or float).")
        if not isinstance(title, str):
            raise ValueError("'title' must be a string.")
        if not isinstance(x_label, str):
            raise ValueError("'x_label' must be a string.")
        if not isinstance(y_label, str):
            raise ValueError("'y_label' must be a string.")
        if num_bins is not None:
            if not isinstance(num_bins, int) or num_bins <= 0:
                raise ValueError("'num_bins' must be a positive integer.")
        if bar_colour is not None and not isinstance(bar_colour, str):
            raise ValueError("'bar_colour' must be a string representing a colour.")
        if not isinstance(margin, int) or margin < 0:
            raise ValueError("'margin' must be a non-negative integer.")
        if not isinstance(filename, str):
            raise ValueError("'filename' must be a string.")
        filename = self._normalize_html_filename(filename)
        if not isinstance(open_browser, bool):
            raise ValueError("'open_browser' must be a Boolean value.")
        
        n = len(variable)
        data_min, data_max = min(variable), max(variable)

        # Calculate bins
        if num_bins is None:
            if n > 30:
                bins = math.ceil(math.log2(n)) + 1  # Sturges' formula
            else:
                bins = len(set(variable))
        else:
            bins = num_bins

        if bins < 2:
            bins = 2
        elif bins > 20:
            bins = 20
        
        num_bins = bins

        # Now use bins to compute bin_edges and counts
        bin_edges = [data_min + i*(data_max-data_min)/bins for i in range(bins+1)]
        counts = [0]*bins
        for value in variable:
            # find bin index
            idx = min(int((value-data_min)/(data_max-data_min)*bins), bins-1)
            counts[idx] += 1

        # Normalize bar colour
        if bar_colour is None:
            hue = random.randint(0, 360)
            bar_colour = f"hsl({hue},70%,50%)"
        
        # Define plot dimensions
        plot_width = 600
        plot_height = 400

        # Define canvas dimensions (plot size + margins)
        canvas_width = 800 + 2*margin
        canvas_height = 600 + 2*margin

        # Begin SVG output
        buf = io.StringIO()
        buf.write(f"<svg xmlns='http://www.w3.org/2000/svg' width='{canvas_width}' height='{canvas_height}' "
                f"viewBox='0 0 {canvas_width} {canvas_height}' style='font-family:sans-serif'>\n")

        # White background rectangle
        buf.write("  <rect width='100%' height='100%' fill='white'/>\n")

        # Title and description for accessibility
        buf.write(f"<title>{title}</title>\n")
        buf.write(f"<desc>Responsive histogram with {num_bins} bins</desc>\n")

        # Plot area (centered)
        plot_left = (canvas_width - plot_width) / 2
        plot_right = plot_left + plot_width
        plot_top = (canvas_height - plot_height) / 2
        plot_bottom = plot_top + plot_height

        # Axes
        buf.write(f"<line x1='{plot_left}' y1='{plot_top}' x2='{plot_left}' y2='{plot_bottom}' stroke='black'/>\n")
        buf.write(f"<line x1='{plot_left}' y1='{plot_bottom}' x2='{plot_right}' y2='{plot_bottom}' stroke='black'/>\n")

        # Title
        buf.write(f"<text x='{(plot_left+plot_right)/2}' y='{plot_top-20}' text-anchor='middle' font-size='16'>{title}</text>\n")

        # Determine max count for scaling
        max_count = max(counts)

        # Coordinate transforms
        def sx(x): return plot_left + (x - data_min) / (data_max - data_min) * (plot_right - plot_left)
        def sy(y): return plot_bottom - (y / max_count) * (plot_bottom - plot_top)

        # Ticks and tick labels
        for i in range(num_bins+1):
            x_pos = sx(bin_edges[i])
            buf.write(f"<line x1='{x_pos}' y1='{plot_bottom}' x2='{x_pos}' y2='{plot_bottom+5}' stroke='black'/>\n")
            buf.write(f"<text x='{x_pos}' y='{plot_bottom+20}' transform='rotate(-45 {x_pos},{plot_bottom+20})' "
                    f"text-anchor='end'>{round(bin_edges[i],2)}</text>\n")

        for j in range(5+1):
            y_val = j*max_count/5
            y_pos = sy(y_val)
            buf.write(f"<line x1='{plot_left-5}' y1='{y_pos}' x2='{plot_left}' y2='{y_pos}' stroke='black'/>\n")
            buf.write(f"<text x='{plot_left-10}' y='{y_pos+4}' text-anchor='end'>{int(y_val)}</text>\n")
        
        # Draw bars
        for i in range(num_bins):
            x0 = sx(bin_edges[i])
            x1 = sx(bin_edges[i+1])
            bar_width = x1 - x0
            y0 = sy(counts[i])
            y1 = plot_bottom
            buf.write(f"<rect x='{x0}' y='{y0}' width='{bar_width}' height='{y1-y0}' "
                    f"fill='{bar_colour}' stroke='black'>\n")
            buf.write(f"<title>Bin {i+1}: {round(bin_edges[i],2)}–{round(bin_edges[i+1],2)}, Count={counts[i]}</title>\n")
            buf.write("</rect>\n")
        
        # Dynamic offset for x-axis label
        max_label_len = max(len(str(round(edge,2))) for edge in bin_edges)
        base_offset = 40
        x_label_offset = base_offset + max_label_len * 5

        # Find longest y tick label
        max_y_label_len = max(len(str(int(j*max_count/5))) for j in range(6))

        # Base offset for small labels
        base_offset = 30

        # Add extra spacing proportional to label length
        y_label_offset = base_offset + max_y_label_len * 6   # tweak multiplier for font size

        # Axis labels (moved here, after ticks)
        buf.write(f"<text x='{(plot_left+plot_right)/2}' y='{plot_bottom+x_label_offset}' "
                f"text-anchor='middle'>{x_label}</text>\n")

        buf.write(f"<text x='{plot_left-y_label_offset}' y='{(plot_top+plot_bottom)/2}' text-anchor='middle' "
            f"transform='rotate(-90 {plot_left-y_label_offset},{(plot_top+plot_bottom)/2})'>{y_label}</text>\n")

        # Close SVG
        buf.write("</svg>\n")

        # Also write HTML wrapper with save button
        html_filename = filename
        html_content = f"""<!DOCTYPE html>
        <html>
        <head>
        <meta charset="utf-8">
        <title>{title}</title>
        <style>
            html, body {{
                height: 100%;
                margin: 0;
            }}
            .wrap {{
                width: 100vw;
                height: 100vh;
                display: flex;
                justify-content: center;
                align-items: center;
            }}
            svg {{
                max-width: 100vw;
                max-height: 100vh;
                width: 100%;
                height: auto;
                display: block;
                background: white;
                font-family: Arial, Helvetica, sans-serif;
            }}
            #saveBtn {{
                position: fixed;
                top: 10px;
                right: 10px;
                padding: 6px 12px;
                font-size: 14px;
                cursor: pointer;
            }}
        </style>
        </head>
        <body>
        <div class="wrap">
            {buf.getvalue()}
        </div>
        <button id="saveBtn" onclick="downloadPNG()">Save as PNG</button>
        <script>
        function downloadPNG() {{
            const svg = document.querySelector("svg");
            const serializer = new XMLSerializer();
            const svgStr = serializer.serializeToString(svg);

            let width, height;
            if (svg.viewBox && svg.viewBox.baseVal) {{
                width = svg.viewBox.baseVal.width;
                height = svg.viewBox.baseVal.height;
            }} else {{
                width = svg.width.baseVal.value;
                height = svg.height.baseVal.value;
            }}

            const canvas = document.createElement("canvas");
            canvas.width = width;
            canvas.height = height;
            const ctx = canvas.getContext("2d");

            const img = new Image();
            img.onload = function() {{
                ctx.drawImage(img, 0, 0, width, height);
                const link = document.createElement("a");
                const safeTitle = "{title}".replace(/[^a-z0-9]+/gi, "_").toLowerCase();
                link.download = safeTitle + ".png";
                link.href = canvas.toDataURL("image/png");
                link.click();
            }};
            img.src = "data:image/svg+xml;base64," + btoa(unescape(encodeURIComponent(svgStr)));
        }}
        </script>
        </body>
        </html>"""

        with open(html_filename, "w", encoding="utf-8") as f:
            f.write(html_content)

        # Optionally open in browser (open HTML wrapper so button is visible)
        if open_browser:
            import webbrowser
            webbrowser.open(html_filename)

    def box_plot(self,data,labels=None,title="Box Plot",x_label="Categories",y_label="Values",num_ticks=5,colour="#4e79a7",margin=60,filename="box_plot.html",open_browser=False):
        """Create a box plot from a list of data series ('data') and configure plot elements. One box will be created per series.

        Parameters:

        data (list): the nested list where each sub-list contains numeric values for one category.

        labels (list): the list of category labels for the x axis. Must be the same length as data.

        title (string): the title of the box plot.

        x_label (string): the label of the x axis.

        y_label (string): the label of the y axis.

        num_ticks (integer): the number of tick marks to draw on the y axis.

        colour (string): the colour of the boxes. If not specified, a random colour will be generated. Colours may be specified as a common colour name or hex code.

        margin (integer): the margin size (in pixels) around the plot area.

        filename (string): the name of the output file. Include the output file path in the string or prepend tt.get_out_dir() to the filename.

        open_browser (Boolean): flag to indicate if the generated file should be automatically opened in the default web browser."""
        if not isinstance(data, list) or not data:
            raise ValueError("'data' must be a non-empty list of data series.")
        if not all(isinstance(series, list) and series for series in data):
            raise ValueError("Each entry in 'data' must be a non-empty list of numeric values.")
        for series in data:
            if not all(isinstance(v, (int, float)) for v in series):
                raise ValueError("All values in each data series must be numeric (int or float).")
        num_series = len(data)
        if labels is not None:
            if not isinstance(labels, list) or len(labels) != num_series:
                raise ValueError("'labels' must be a list with one label per data series.")
            if not all(isinstance(lbl, str) for lbl in labels):
                raise ValueError("All entries in 'labels' must be strings.")
        if not isinstance(title, str):
            raise ValueError("'title' must be a string.")
        if not isinstance(x_label, str):
            raise ValueError("'x_label' must be a string.")
        if not isinstance(y_label, str):
            raise ValueError("'y_label' must be a string.")
        if not isinstance(num_ticks, int) or num_ticks <= 0:
            raise ValueError("'num_ticks' must be a positive integer.")
        if colour is not None and not isinstance(colour, str):
            raise ValueError("'colour' must be a string representing a colour.")
        if not isinstance(margin, int) or margin < 0:
            raise ValueError("'margin' must be a non-negative integer.")
        if not isinstance(filename, str):
            raise ValueError("'filename' must be a string.")
        filename = self._normalize_html_filename(filename)
        if not isinstance(open_browser, bool):
            raise ValueError("'open_browser' must be a Boolean value.")

        if labels is None:
            labels = [f"Category {i+1}" for i in range(len(data))]
        else:
            labels = [str(lbl) for lbl in labels]  # ensure string labels

        n_cat = len(data)

        # Dynamic margins (same pattern as bar chart)
        extra_left = margin + 40   # space for y-axis labels
        extra_right = margin
        extra_top = margin + 40    # space for title
        extra_bottom = margin + 60 # space for x-axis labels

        width = 800 + extra_left + extra_right
        height = 600 + extra_top + extra_bottom

        plot_left = extra_left
        plot_right = width - extra_right
        plot_top = extra_top
        plot_bottom = height - extra_bottom

        # Axis ranges
        all_values = [v for series in data for v in series]
        y_min = min(all_values)
        y_max = max(all_values)

        # Slot width: total plot width divided by number of categories
        slot_width = (plot_right - plot_left) / n_cat

        # Optimal box width: leave some spacing between boxes
        box_width = slot_width * 0.6   # 60% of slot width, tweak factor as needed

        # Quartile and whisker calculations
        stats = []
        for series in data:
            s = sorted(series)
            # Quartiles (Q1, Median, Q3)
            if len(s) >= 2:
                q = statistics.quantiles(s, n=4, method="exclusive")
                q1, q3 = q[0], q[2]
            else:
                q1, q3 = s[0], s[0]
            median = statistics.median(s)
            iqr = q3 - q1

            # Whiskers (Tukey rule: 1.5 * IQR)
            lower_cut = q1 - 1.5 * iqr
            upper_cut = q3 + 1.5 * iqr
            lower_whisker = min((v for v in s if v >= lower_cut), default=q1)
            upper_whisker = max((v for v in s if v <= upper_cut), default=q3)

            # Outliers
            outliers = [v for v in s if v < lower_whisker or v > upper_whisker]

            stats.append({
                "q1": q1,
                "median": median,
                "q3": q3,
                "lower_whisker": lower_whisker,
                "upper_whisker": upper_whisker,
                "outliers": outliers
            })

        # Scaling function (map values to SVG coordinates)
        def scale_y(val):
            return plot_bottom - (val - y_min) / (y_max - y_min) * (plot_bottom - plot_top)

            # === Begin SVG output ===
        buffer = io.StringIO()
        buffer.write(
            f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' "
            f"viewBox='0 0 {width} {height}' style='font-family:sans-serif'>\n"
        )
        buffer.write('  <rect width="100%" height="100%" fill="white"/>\n')

        # Title and description for accessibility
        buffer.write(f"<title>{title}</title>\n")
        buffer.write(f"<desc>Box plot with {n_cat} categories</desc>\n")

        # Chart title
        buffer.write(
            f"<text x='{(plot_left+plot_right)/2}' y='{margin}' "
            f"text-anchor='middle' font-size='16'>{title}</text>\n"
        )

        # Y axis line
        buffer.write(f"<line x1='{plot_left}' y1='{plot_top}' x2='{plot_left}' y2='{plot_bottom}' stroke='black'/>\n")
        # X axis line
        buffer.write(f"<line x1='{plot_left}' y1='{plot_bottom}' x2='{plot_right}' y2='{plot_bottom}' stroke='black'/>\n")

        # Y axis ticks and labels
        for i in range(num_ticks + 1):
            val = y_min + i * (y_max - y_min) / num_ticks
            y = scale_y(val)
            buffer.write(f"<line x1='{plot_left-5}' y1='{y}' x2='{plot_left}' y2='{y}' stroke='black'/>\n")
            buffer.write(f"<text x='{plot_left-10}' y='{y+4}' text-anchor='end' font-size='12'>{round(val,2)}</text>\n")

        # Dynamic offsets for axis labels
        max_label_len = max(len(str(lbl)) for lbl in labels)
        x_label_offset = 40 + max_label_len * 5
        max_y_label_len = max(len(str(round(y_min + i * (y_max - y_min) / num_ticks,2))) for i in range(num_ticks+1))
        y_label_offset = 30 + max_y_label_len * 6

        # Axis labels
        buffer.write(
            f"<text x='{(plot_left+plot_right)/2}' y='{plot_bottom+x_label_offset}' "
            f"text-anchor='middle' font-size='14'>{x_label}</text>\n"
        )
        buffer.write(
            f"<text x='{plot_left-y_label_offset}' y='{(plot_top+plot_bottom)/2}' text-anchor='middle' "
            f"font-size='14' transform='rotate(-90 {plot_left-y_label_offset},{(plot_top+plot_bottom)/2})'>{y_label}</text>\n"
        )

        # Draw boxes, whiskers, medians, outliers
        for i, stat in enumerate(stats):
            x_center = plot_left + (i + 0.5) * slot_width
            half_width = box_width / 2

            y_q1 = scale_y(stat["q1"])
            y_q3 = scale_y(stat["q3"])
            y_median = scale_y(stat["median"])
            y_low = scale_y(stat["lower_whisker"])
            y_high = scale_y(stat["upper_whisker"])

            # Box (Q1–Q3)
            buffer.write(
                f"<rect x='{x_center-half_width}' y='{y_q3}' width='{box_width}' height='{y_q1-y_q3}' "
                f"fill='{colour}' stroke='black'>\n"
            )
            buffer.write(f"<title>{labels[i]} Q1={stat['q1']}, Median={stat['median']}, Q3={stat['q3']}</title>\n")
            buffer.write("</rect>\n")

            # Median line
            buffer.write(
                f"<line x1='{x_center-half_width}' y1='{y_median}' x2='{x_center+half_width}' y2='{y_median}' "
                f"stroke='black' stroke-width='2'><title>{labels[i]} Median={stat['median']}</title></line>\n"
            )

            # Whiskers
            buffer.write(f"<line x1='{x_center}' y1='{y_q3}' x2='{x_center}' y2='{y_high}' stroke='black'/>\n")
            buffer.write(f"<line x1='{x_center}' y1='{y_q1}' x2='{x_center}' y2='{y_low}' stroke='black'/>\n")
            buffer.write(f"<line x1='{x_center-half_width/2}' y1='{y_high}' x2='{x_center+half_width/2}' y2='{y_high}' stroke='black'/>\n")
            buffer.write(f"<line x1='{x_center-half_width/2}' y1='{y_low}' x2='{x_center+half_width/2}' y2='{y_low}' stroke='black'/>\n")

            # Outliers
            for out in stat["outliers"]:
                y_out = scale_y(out)
                buffer.write(f"<circle cx='{x_center}' cy='{y_out}' r='3' fill='black'><title>{labels[i]} Outlier={out}</title></circle>\n")

            # X axis category labels
            buffer.write(
                f"<text x='{x_center}' y='{plot_bottom+20}' text-anchor='end' font-size='12' "
                f"transform='rotate(-45 {x_center},{plot_bottom+20})'>{labels[i]}</text>\n"
            )

        # Close SVG
        buffer.write("</svg>\n")

        # Write HTML wrapper with save button
        html_filename = filename
        html_content = f"""<!DOCTYPE html>
        <html>
        <head>
        <meta charset="utf-8">
        <title>{title}</title>
        <style>
            html, body {{
                height: 100%;
                margin: 0;
            }}
            .wrap {{
                width: 100vw;
                height: 100vh;
                display: flex;
                justify-content: center;
                align-items: center;
            }}
            svg {{
                max-width: 100vw;
                max-height: 100vh;
                width: 100%;
                height: auto;
                display: block;
                background: white;
                font-family: Arial, Helvetica, sans-serif;
            }}
            #saveBtn {{
                position: fixed;
                top: 10px;
                right: 10px;
                padding: 6px 12px;
                font-size: 14px;
                cursor: pointer;
            }}
        </style>
        </head>
        <body>
        <div class="wrap">
            {buffer.getvalue()}
        </div>
        <button id="saveBtn" onclick="downloadPNG()">Save as PNG</button>
        <script>
        function downloadPNG() {{
            const svg = document.querySelector("svg");
            const serializer = new XMLSerializer();
            const svgStr = serializer.serializeToString(svg);

            let width, height;
            if (svg.viewBox && svg.viewBox.baseVal) {{
                width = svg.viewBox.baseVal.width;
                height = svg.viewBox.baseVal.height;
            }} else {{
                width = svg.width.baseVal.value;
                height = svg.height.baseVal.value;
            }}

            const canvas = document.createElement("canvas");
            canvas.width = width;
            canvas.height = height;
            const ctx = canvas.getContext("2d");

            const img = new Image();
            img.onload = function() {{
                ctx.drawImage(img, 0, 0, width, height);
                const link = document.createElement("a");
                const safeTitle = "{title}".replace(/[^a-z0-9]+/gi, "_").toLowerCase();
                link.download = safeTitle + ".png";
                link.href = canvas.toDataURL("image/png");
                link.click();
            }};
            img.src = "data:image/svg+xml;base64," + btoa(unescape(encodeURIComponent(svgStr)));
        }}
        </script>
        </body>
        </html>"""

        with open(html_filename, "w", encoding="utf-8") as f:
            f.write(html_content)

        # Optionally open in browser
        if open_browser:
            import webbrowser
            webbrowser.open(html_filename)

    def grouped_box_plot(self,data,group_labels,category_labels,title="Grouped Box Plot",x_label="Groups",y_label="Values",num_ticks=5,colours=None,margin=60,filename="grouped_box_plot.html",open_browser=False):
        """Create a grouped box plot for multiple categories across groups and configure plot elements.

        Parameters:

        data (list): the nested list representing groups and categories to visualize. The outermost lists are groups, the mid-level lists are categories within groups, and the innermost lists must contain numeric values.
                            
        group_labels (list): labels for each group (e.g., years). Must match the number of groups in data.

        category_labels (list): labels for categories within each group. Must match the number of categories in each group.

        title (string): the title of the box plot.

        x_label (string): the label of the x axis.

        y_label (string): the label of the y axis.

        num_ticks (integer): the number of tick marks to draw on the y axis.

        colours (list): the colour of the boxes, one per category. If not specified, a random colour will be generated. Colours may be specified as a common colour name or hex code.

        margin (integer): the margin size (in pixels) around the plot area.

        filename (string): the name of the output file. Include the output file path in the string or prepend tt.get_out_dir() to the filename.

        open_browser (Boolean): flag to indicate if the generated file should be automatically opened in the default web browser."""
        if not isinstance(data, list) or not data:
            raise ValueError("'data' must be a non-empty list of groups.")
        if not all(isinstance(group, list) and group for group in data):
            raise ValueError("Each group in 'data' must be a non-empty list of categories.")
        num_groups = len(data)
        num_categories = None
        for group in data:
            if not all(isinstance(category, list) and category for category in group):
                raise ValueError("Each category in each group must be a non-empty list of numeric values.")
            if not all(all(isinstance(v, (int, float)) for v in category) for category in group):
                raise ValueError("All values in each category must be numeric (int or float).")
            if num_categories is None:
                num_categories = len(group)
            elif len(group) != num_categories:
                raise ValueError("All groups must contain the same number of categories.")
        if not isinstance(group_labels, list) or len(group_labels) != num_groups:
            raise ValueError("'group_labels' must be a list with one label per group.")
        if not all(isinstance(lbl, str) for lbl in group_labels):
            raise ValueError("All entries in 'group_labels' must be strings.")
        if not isinstance(category_labels, list) or len(category_labels) != num_categories:
            raise ValueError("'category_labels' must be a list with one label per category.")
        if not all(isinstance(lbl, str) for lbl in category_labels):
            raise ValueError("All entries in 'category_labels' must be strings.")
        if not isinstance(title, str):
            raise ValueError("'title' must be a string.")
        if not isinstance(x_label, str):
            raise ValueError("'x_label' must be a string.")
        if not isinstance(y_label, str):
            raise ValueError("'y_label' must be a string.")
        if not isinstance(num_ticks, int) or num_ticks <= 0:
            raise ValueError("'num_ticks' must be a positive integer.")
        if colours is not None:
            if not isinstance(colours, list):
                raise ValueError("'colours' must be a list of colour strings.")
            if len(colours) != num_categories:
                raise ValueError("Length of 'colours' must match the number of categories.")
            if not all(isinstance(c, str) for c in colours):
                raise ValueError("All entries in 'colours' must be strings representing colours.")
        if not isinstance(margin, int) or margin < 0:
            raise ValueError("'margin' must be a non-negative integer.")
        if not isinstance(filename, str):
            raise ValueError("'filename' must be a string.")
        filename = self._normalize_html_filename(filename)
        if not isinstance(open_browser, bool):
            raise ValueError("'open_browser' must be a Boolean value.")        

        n_groups = len(data)
        n_categories = len(category_labels)

        # Assign colours
        if colours is None:
            base_step = 360 / n_categories
            offset = random.randint(0, 359)  # random starting hue
            colours = [
                f"hsl({int((offset + i * base_step) % 360)},70%,50%)"
                for i in range(n_categories)
            ]
        elif len(colours) < n_categories:
            raise ValueError("Not enough colours provided for categories.")
        
        # Fixed canvas dimensions
        width = 800
        height = 600
        margin_left = margin * 2
        margin_top = margin

        # Axis ranges (based on all numeric values across groups/categories)
        all_values = [v for group in data for series in group for v in series]
        y_min = min(all_values)
        y_max = max(all_values)

        # Dynamic axis offsets based on longest tick labels
        longest_x_len = max(len(str(gl)) for gl in group_labels) if n_groups > 0 else 0
        longest_y_len = max(
            len(str(round(y_min + i * (y_max - y_min) / num_ticks, 2)))
            for i in range(num_ticks + 1)
        )
        x_label_offset = 40 + longest_x_len * 5
        y_label_offset = 30 + longest_y_len * 6

        # Legend sizing
        legend_font = 12
        swatch_size = 14
        swatch_gap = 10
        legend_row_gap = 6
        est_text_widths = [len(str(cat)) * 8 for cat in category_labels]
        legend_inner_width = max(swatch_size + swatch_gap + tw for tw in est_text_widths)

        # Ensure the legend title always fits — set a minimum width
        min_title_width = len("Legend") * 8 + 20  # rough estimate: 8px per char + padding
        legend_width = max(legend_inner_width + 20, min_title_width)

        legend_height = n_categories * (swatch_size + legend_row_gap) + 30



        # Plot area (centered, smaller footprint to leave room for labels/legend)
        plot_width = width * 0.65
        plot_height = height * 0.65
        plot_left = (width - plot_width) / 2
        plot_right = plot_left + plot_width
        plot_top = (height - plot_height) / 2
        plot_bottom = plot_top + plot_height

        # Group slot geometry with automatic box sizing
        slot_width = (plot_right - plot_left) / n_groups
        group_gap_fraction = 0.2  # reserve 20% of slot for separation
        box_area_width = slot_width * (1 - group_gap_fraction)
        box_width_px = box_area_width / n_categories
        total_box_block = n_categories * box_width_px

        # Legend position (outside plot, inside canvas)
        legend_x = plot_right + 40
        legend_y = plot_top

        buffer = io.StringIO()
        buffer.write(
            f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' "
            f"viewBox='0 0 {width} {height}' style='font-family:sans-serif'>\n"
        )
        buffer.write("  <rect width='100%' height='100%' fill='white'/>\n")

        # Accessibility
        buffer.write(f"<title>{title}</title>\n")
        buffer.write(f"<desc>Grouped box plot with {n_groups} groups and {n_categories} categories</desc>\n")

        # Title
        buffer.write(
            f"<text x='{(plot_left+plot_right)/2}' y='{plot_top-20}' "
            f"text-anchor='middle' font-size='16'>{title}</text>\n"
        )

        # Axes
        buffer.write(f"<line x1='{plot_left}' y1='{plot_top}' x2='{plot_left}' y2='{plot_bottom}' stroke='black'/>\n")
        buffer.write(f"<line x1='{plot_left}' y1='{plot_bottom}' x2='{plot_right}' y2='{plot_bottom}' stroke='black'/>\n")

        # Y ticks
        def scale_y(val):
            return plot_bottom - (val - y_min) / (y_max - y_min) * (plot_bottom - plot_top)

        for i in range(num_ticks + 1):
            val = y_min + i * (y_max - y_min) / num_ticks
            y = scale_y(val)
            buffer.write(f"<line x1='{plot_left-5}' y1='{y}' x2='{plot_left}' y2='{y}' stroke='black'/>\n")
            buffer.write(f"<text x='{plot_left-10}' y='{y+4}' text-anchor='end' font-size='12'>{round(val,2)}</text>\n")

        # Axis labels
        buffer.write(
            f"<text x='{(plot_left+plot_right)/2}' y='{plot_bottom + x_label_offset}' "
            f"text-anchor='middle' font-size='14'>{x_label}</text>\n"
        )
        buffer.write(
            f"<text x='{plot_left - y_label_offset}' y='{(plot_top+plot_bottom)/2}' text-anchor='middle' "
            f"font-size='14' transform='rotate(-90 {plot_left - y_label_offset},{(plot_top+plot_bottom)/2})'>{y_label}</text>\n"
        )

        # Group labels
        for g, glabel in enumerate(group_labels):
            group_center = plot_left + (g + 0.5) * slot_width
            buffer.write(
                f"<text x='{group_center}' y='{plot_bottom + 20}' font-size='12' text-anchor='end' "
                f"transform='rotate(-45 {group_center},{plot_bottom + 20})'>{glabel}</text>\n"
            )
        
        # Rendering grouped boxes
        for g, group in enumerate(data):
            group_left = plot_left + g * slot_width
            start_x = group_left + (slot_width - box_area_width) / 2

            for c, series in enumerate(group):
                s = sorted(series)
                # Quartiles
                if len(s) >= 2:
                    q = statistics.quantiles(s, n=4, method="exclusive")
                    q1, q3 = q[0], q[2]
                else:
                    q1, q3 = s[0], s[0]
                median = statistics.median(s)
                iqr = q3 - q1

                # Whiskers
                lower_cut = q1 - 1.5 * iqr
                upper_cut = q3 + 1.5 * iqr
                lower_whisker = min((v for v in s if v >= lower_cut), default=q1)
                upper_whisker = max((v for v in s if v <= upper_cut), default=q3)

                # Outliers
                outliers = [v for v in s if v < lower_whisker or v > upper_whisker]

                # Coordinates
                x_left = start_x + c * box_width_px
                x_center = x_left + box_width_px / 2
                half_width = box_width_px / 2

                y_q1 = scale_y(q1)
                y_q3 = scale_y(q3)
                y_median = scale_y(median)
                y_low = scale_y(lower_whisker)
                y_high = scale_y(upper_whisker)

                # Box (Q1–Q3)
                buffer.write(
                    f"<rect x='{x_center-half_width}' y='{y_q3}' width='{box_width_px}' height='{y_q1-y_q3}' "
                    f"fill='{colours[c]}' stroke='black'>\n"
                )
                buffer.write(f"<title>{group_labels[g]} – {category_labels[c]} – Q1={q1}, Median={median}, Q3={q3}</title>\n")
                buffer.write("</rect>\n")

                # Median line
                buffer.write(
                    f"<line x1='{x_center-half_width}' y1='{y_median}' x2='{x_center+half_width}' y2='{y_median}' "
                    f"stroke='black' stroke-width='2'><title>{group_labels[g]} – {category_labels[c]} Median={median}</title></line>\n"
                )

                # Whiskers
                buffer.write(f"<line x1='{x_center}' y1='{y_q3}' x2='{x_center}' y2='{y_high}' stroke='black'/>\n")
                buffer.write(f"<line x1='{x_center}' y1='{y_q1}' x2='{x_center}' y2='{y_low}' stroke='black'/>\n")
                buffer.write(f"<line x1='{x_center-half_width/2}' y1='{y_high}' x2='{x_center+half_width/2}' y2='{y_high}' stroke='black'/>\n")
                buffer.write(f"<line x1='{x_center-half_width/2}' y1='{y_low}' x2='{x_center+half_width/2}' y2='{y_low}' stroke='black'/>\n")

                # Outliers
                for out in outliers:
                    y_out = scale_y(out)
                    buffer.write(f"<circle cx='{x_center}' cy='{y_out}' r='3' fill='black'><title>{group_labels[g]} – {category_labels[c]} Outlier={out}</title></circle>\n")

        # Legend
        buffer.write(
            f"<rect x='{legend_x}' y='{legend_y}' width='{legend_width}' height='{legend_height}' "
            f"fill='white' stroke='black' rx='6' ry='6'/>\n"
        )
        buffer.write(
            f"<text x='{legend_x + 10}' y='{legend_y + 18}' font-size='14' text-anchor='start'>Legend</text>\n"
        )
        for i, cat in enumerate(category_labels):
            entry_y = legend_y + 30 + legend_row_gap + i * (swatch_size + legend_row_gap)
            buffer.write(
                f"<rect x='{legend_x + 10}' y='{entry_y - swatch_size + 2}' "
                f"width='{swatch_size}' height='{swatch_size}' fill='{colours[i]}' stroke='black'/>\n"
            )
            buffer.write(
                f"<text x='{legend_x + 10 + swatch_size + swatch_gap}' y='{entry_y}' "
                f"font-size='{legend_font}' text-anchor='start'>{cat}</text>\n"
            )

        buffer.write("</svg>\n")

        # HTML wrapper with save button
        html_filename = filename
        html_content = f"""<!DOCTYPE html>
        <html>
        <head>
        <meta charset="utf-8">
        <title>{title}</title>
        <style>
            html, body {{
                height: 100%;
                margin: 0;
            }}
            .wrap {{
                width: 100vw;
                height: 100vh;
                display: flex;
                justify-content: center;
                align-items: center;
            }}
            svg {{
                max-width: 100vw;
                max-height: 100vh;
                width: 100%;
                height: auto;
                display: block;
                background: white;
                font-family: Arial, Helvetica, sans-serif;
            }}
            #saveBtn {{
                position: fixed;
                top: 10px;
                right: 10px;
                padding: 6px 12px;
                font-size: 14px;
                cursor: pointer;
            }}
        </style>
        </head>
        <body>
        <div class="wrap">
            {buffer.getvalue()}
        </div>
        <button id="saveBtn" onclick="downloadPNG()">Save as PNG</button>
        <script>
        function downloadPNG() {{
            const svg = document.querySelector("svg");
            const serializer = new XMLSerializer();
            const svgStr = serializer.serializeToString(svg);

            let width, height;
            if (svg.viewBox && svg.viewBox.baseVal) {{
                width = svg.viewBox.baseVal.width;
                height = svg.viewBox.baseVal.height;
            }} else {{
                width = svg.width.baseVal.value;
                height = svg.height.baseVal.value;
            }}

            const canvas = document.createElement("canvas");
            canvas.width = width;
            canvas.height = height;
            const ctx = canvas.getContext("2d");

            const img = new Image();
            img.onload = function() {{
                ctx.drawImage(img, 0, 0, width, height);
                const link = document.createElement("a");
                const safeTitle = "{title}".replace(/[^a-z0-9]+/gi, "_").toLowerCase();
                link.download = safeTitle + ".png";
                link.href = canvas.toDataURL("image/png");
                link.click();
            }};
            img.src = "data:image/svg+xml;base64," + btoa(unescape(encodeURIComponent(svgStr)));
        }}
        </script>
        </body>
        </html>"""

        with open(html_filename, "w", encoding="utf-8") as f:
            f.write(html_content)

        if open_browser:
            import webbrowser
            webbrowser.open(html_filename)

    # ============================================================
    # CATEGORICAL PLOTS
    # ============================================================

    def bar_chart(self,variable,title="Bar Chart",x_label="Category",y_label="Count",num_ticks=6,colour="#4e79a7",margin=60,filename="bar_chart.html",open_browser=False):
        """Create a bar chart from a list of categorical values ('variable') and configure chart elements.

        Parameters:

        variable (list): the list of categorical values. 

        title (string): the title of the bar chart.

        x_label (string): the label of the x axis.

        y_label (string): the label of the y axis.

        num_ticks (integer): The number of tick marks to draw on the y axis.

        colour (string): the colour of the bars. If not specified, a random colour will be generated. Colours may be specified as a common colour name or hex code.

        margin (integer): the margin size (in pixels) around the plot area.

        filename (string): the name of the output file. Include the output file path in the string or prepend tt.get_out_dir() to the filename.

        open_browser (Boolean): flag to indicate if the generated file should be automatically opened in the default web browser."""
        if not isinstance(variable, list) or not variable:
            raise ValueError("'variable' must be a non-empty list of categorical values.")
        if not isinstance(title, str):
            raise ValueError("'title' must be a string.")
        if not isinstance(x_label, str):
            raise ValueError("'x_label' must be a string.")
        if not isinstance(y_label, str):
            raise ValueError("'y_label' must be a string.")
        if not isinstance(num_ticks, int) or num_ticks <= 0:
            raise ValueError("'num_ticks' must be a positive integer.")
        if colour is not None and not isinstance(colour, str):
            raise ValueError("'colour' must be a string representing a colour.")
        if not isinstance(margin, int) or margin < 0:
            raise ValueError("'margin' must be a non-negative integer.")
        if not isinstance(filename, str):
            raise ValueError("'filename' must be a string.")
        filename = self._normalize_html_filename(filename)
        if not isinstance(open_browser, bool):
            raise ValueError("'open_browser' must be a Boolean value.")

        from collections import Counter
        # Count categories
        counts = Counter(variable)
        categories = list(counts.keys())
        n_cat = len(categories)
        total = sum(counts.values())

        # Dynamic margins (similar to line/histogram plots)
        extra_left = margin + 40   # space for y-axis labels
        extra_right = margin
        extra_top = margin + 40    # space for title
        extra_bottom = margin + 60 # space for x-axis labels

        width = 800 + extra_left + extra_right
        height = 600 + extra_top + extra_bottom

        plot_left = extra_left
        plot_right = width - extra_right
        plot_top = extra_top
        plot_bottom = height - extra_bottom

        # Axis ranges
        max_count = max(counts.values())
        y_max = max_count
        y_min = 0

        # Bar positioning (adaptive width)
        slot_width = (plot_right - plot_left) / n_cat
        bar_width = slot_width * 0.7   # bars take 70% of slot
        gap = slot_width - bar_width   # remaining space is gap

        # Begin SVG output
        buffer = io.StringIO()
        buffer.write(
            f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' "
            f"viewBox='0 0 {width} {height}' style='font-family:sans-serif'>\n"
        )
        buffer.write('  <rect width="100%" height="100%" fill="white"/>\n')

        # Title and description for accessibility
        buffer.write(f"<title>{title}</title>\n")
        buffer.write(f"<desc>Bar chart with {n_cat} categories</desc>\n")

        # Chart title
        buffer.write(
            f"<text x='{(plot_left+plot_right)/2}' y='{margin}' "
            f"text-anchor='middle' font-size='16'>{title}</text>\n"
        )

        # Y axis line
        buffer.write(f"<line x1='{plot_left}' y1='{plot_top}' x2='{plot_left}' y2='{plot_bottom}' stroke='black'/>\n")
        # X axis line
        buffer.write(f"<line x1='{plot_left}' y1='{plot_bottom}' x2='{plot_right}' y2='{plot_bottom}' stroke='black'/>\n")

        # Y axis ticks and labels
        for i in range(num_ticks + 1):
            val = y_min + i * (y_max - y_min) / num_ticks
            y = plot_bottom - (val - y_min) / (y_max - y_min) * (plot_bottom - plot_top)
            buffer.write(f"<line x1='{plot_left-5}' y1='{y}' x2='{plot_left}' y2='{y}' stroke='black'/>\n")
            buffer.write(f"<text x='{plot_left-10}' y='{y+4}' text-anchor='end' font-size='12'>{int(val)}</text>\n")

        # Dynamic offset for x-axis label
        max_label_len = max(len(str(cat)) for cat in categories)
        base_offset = 40
        x_label_offset = base_offset + max_label_len * 5

        # Find longest y tick label
        max_y_label_len = max(len(str(int(y_min + i * (y_max - y_min) / num_ticks))) for i in range(num_ticks+1))
        base_offset = 30
        y_label_offset = base_offset + max_y_label_len * 6   # tweak multiplier for font size

        # Axis labels
        buffer.write(
            f"<text x='{(plot_left+plot_right)/2}' y='{plot_bottom+x_label_offset}' "
            f"text-anchor='middle' font-size='14'>{x_label}</text>\n"
        )
        buffer.write(
            f"<text x='{plot_left-y_label_offset}' y='{(plot_top+plot_bottom)/2}' text-anchor='middle' "
            f"font-size='14' transform='rotate(-90 {plot_left-y_label_offset},{(plot_top+plot_bottom)/2})'>{y_label}</text>\n"
        )

        # Draw bars
        for i, category in enumerate(categories):
            count = counts[category]
            frac = count / total
            bar_height = (count - y_min) / (y_max - y_min) * (plot_bottom - plot_top)
            x_center = plot_left + (i + 0.5) * slot_width
            x_left = x_center - bar_width / 2
            y_top = plot_bottom - bar_height

            buffer.write(
                f"<rect x='{x_left}' y='{y_top}' width='{bar_width}' height='{bar_height}' "
                f"fill='{colour}' stroke='black'>\n"
            )
            # Tooltip
            buffer.write(
                f"  <title>{category} – {count} ({frac:.1%})</title>\n"
            )
            buffer.write("</rect>\n")

            # X axis category labels
            buffer.write(
                f"<text x='{x_center}' y='{plot_bottom+20}' text-anchor='end' font-size='12' "
                f"transform='rotate(-45 {x_center},{plot_bottom+20})'>{category}</text>\n"
            )

        # Close SVG
        buffer.write("</svg>\n")

        # Write HTML wrapper with save button
        html_filename = filename
        html_content = f"""<!DOCTYPE html>
        <html>
        <head>
        <meta charset="utf-8">
        <title>{title}</title>
        <style>
            html, body {{
                height: 100%;
                margin: 0;
            }}
            .wrap {{
                width: 100vw;
                height: 100vh;
                display: flex;
                justify-content: center;
                align-items: center;
            }}
            svg {{
                max-width: 100vw;
                max-height: 100vh;
                width: 100%;
                height: auto;
                display: block;
                background: white;
                font-family: Arial, Helvetica, sans-serif;
            }}
            #saveBtn {{
                position: fixed;
                top: 10px;
                right: 10px;
                padding: 6px 12px;
                font-size: 14px;
                cursor: pointer;
            }}
        </style>
        </head>
        <body>
        <div class="wrap">
            {buffer.getvalue()}
        </div>
        <button id="saveBtn" onclick="downloadPNG()">Save as PNG</button>
        <script>
        function downloadPNG() {{
            const svg = document.querySelector("svg");
            const serializer = new XMLSerializer();
            const svgStr = serializer.serializeToString(svg);

            let width, height;
            if (svg.viewBox && svg.viewBox.baseVal) {{
                width = svg.viewBox.baseVal.width;
                height = svg.viewBox.baseVal.height;
            }} else {{
                width = svg.width.baseVal.value;
                height = svg.height.baseVal.value;
            }}

            const canvas = document.createElement("canvas");
            canvas.width = width;
            canvas.height = height;
            const ctx = canvas.getContext("2d");

            const img = new Image();
            img.onload = function() {{
                ctx.drawImage(img, 0, 0, width, height);
                const link = document.createElement("a");
                const safeTitle = "{title}".replace(/[^a-z0-9]+/gi, "_").toLowerCase();
                link.download = safeTitle + ".png";
                link.href = canvas.toDataURL("image/png");
                link.click();
            }};
            img.src = "data:image/svg+xml;base64," + btoa(unescape(encodeURIComponent(svgStr)));
        }}
        </script>
        </body>
        </html>"""

        with open(html_filename, "w", encoding="utf-8") as f:
            f.write(html_content)

        # Optionally open in browser
        if open_browser:
            import webbrowser
            webbrowser.open(html_filename)
   
    def grouped_bar_chart(self,variable,group_labels,category_labels=None,title="Grouped Bar Chart",x_label="Groups",y_label="Count",num_ticks=2,colours=None,margin=60,filename="grouped_bar_chart.html",open_browser=False):
        """Create a grouped bar chart for categorical or numerical input data and configure plot elements.
            
        Parameters:
            
        variable (list): the list of values representing the data to be visualized. If values are strings, categories will be inferred automatically. If values are numerical counts, then category_labels must be supplied.
            
        group_labels (list): the labels for each group. Must match the number of sub-lists in variable.
            
        category_labels (list): the labels for categories within each group. Required if values are numerical counts. Ignored if values are categorical strings, in which case categories are inferred automatically.
            
        title (string): the title of the bar chart.
            
        x_label (string): the label of the x axis.
            
        y_label (string): the label of the y axis.
            
        num_ticks (integer): the number of tick marks to draw on the y axis.
            
        colours (list): the colour of the bars, one per category. If not specified, a random colour will be generated. Colours may be specified as a common colour name or hex code.
            
        margin (integer): the margin size (in pixels) around the plot area.

        filename (string): the name of the output file. Include the output file path in the string or prepend tt.get_out_dir() to the filename.

        open_browser (Boolean): flag to indicate if the generated file should be automatically opened in the default web browser."""
        if not isinstance(variable, list) or not variable:
            raise ValueError("'variable' must be a non-empty list of groups.")
        if not all(isinstance(group, list) and group for group in variable):
            raise ValueError("Each entry in 'variable' must be a non-empty list representing a group.")
        is_all_numeric = all(all(isinstance(v, (int, float)) for v in group) for group in variable)
        is_all_string = all(all(isinstance(v, str) for v in group) for group in variable)
        if not is_all_numeric and not is_all_string:
            raise ValueError("'variable' must contain either all numeric values or all string values within each group.")
        num_groups = len(variable)
        num_categories = len(variable[0])
        if not isinstance(group_labels, list) or len(group_labels) != num_groups:
            raise ValueError("'group_labels' must be a list with one label per group.")
        if not all(isinstance(lbl, str) for lbl in group_labels):
            raise ValueError("All entries in 'group_labels' must be strings.")
        if is_all_numeric:
            if category_labels is None:
                raise ValueError("'category_labels' must be provided when 'variable' contains numeric counts.")
            if not isinstance(category_labels, list) or len(category_labels) != num_categories:
                raise ValueError("'category_labels' must be a list with one label per category.")
            if not all(isinstance(lbl, str) for lbl in category_labels):
                raise ValueError("All entries in 'category_labels' must be strings.")
        else:
            if category_labels is not None:
                if not isinstance(category_labels, list) or not all(isinstance(lbl, str) for lbl in category_labels):
                    raise ValueError("'category_labels' must be a list of strings if provided.")
        if not isinstance(title, str):
            raise ValueError("'title' must be a string.")
        if not isinstance(x_label, str):
            raise ValueError("'x_label' must be a string.")
        if not isinstance(y_label, str):
            raise ValueError("'y_label' must be a string.")
        if not isinstance(num_ticks, int) or num_ticks <= 0:
            raise ValueError("'num_ticks' must be a positive integer.")
        if colours is not None:
            if not isinstance(colours, list):
                raise ValueError("'colours' must be a list of colour strings.")
            if len(colours) != num_categories:
                raise ValueError("Length of 'colours' must match the number of categories.")
            if not all(isinstance(c, str) for c in colours):
                raise ValueError("All entries in 'colours' must be strings representing colours.")
        if not isinstance(margin, int) or margin < 0:
            raise ValueError("'margin' must be a non-negative integer.")
        if not isinstance(filename, str):
            raise ValueError("'filename' must be a string.")
        filename = self._normalize_html_filename(filename)
        if not isinstance(open_browser, bool):
            raise ValueError("'open_browser' must be a Boolean value.")

        # Detect input type: categorical strings vs numeric counts
        first_elem = variable[0][0]
        if isinstance(first_elem, str):
            # Case A: raw categorical values
            categories = sorted(set(cat for group in variable for cat in group))
            counts = []
            for group in variable:
                group_count = {cat: group.count(cat) for cat in categories}
                counts.append(group_count)
        elif isinstance(first_elem, (int, float)):
            # Case B: numeric counts
            if category_labels is None:
                raise ValueError("Numeric input requires category_labels.")
            categories = list(category_labels)
            for idx, group in enumerate(variable):
                if len(group) != len(categories):
                    raise ValueError(
                        f"Group {idx} length ({len(group)}) must match len(category_labels) ({len(categories)})."
                    )
            counts = [{cat: val for cat, val in zip(categories, group)} for group in variable]
        else:
            raise ValueError("Values must be either categorical strings or numeric counts.")

        n_groups = len(variable)
        n_categories = len(categories)

        # Assign colours
        if colours is None:
            base_step = 360 / n_categories
            offset = random.randint(0, 359)  # random starting hue
            colours = [
                f"hsl({int((offset + i * base_step) % 360)},70%,50%)"
                for i in range(n_categories)
            ]
        elif len(colours) < n_categories:
            raise ValueError("Not enough colours provided for categories.")

        # Fixed canvas dimensions
        width = 800
        height = 600
        margin_left = margin * 2
        margin_top = margin

        # Axis ranges
        y_min = 0
        y_max = max(max(group_count.values()) for group_count in counts) or 1

        # Dynamic axis offsets based on longest tick labels
        longest_x_len = max(len(str(gl)) for gl in group_labels) if n_groups > 0 else 0
        longest_y_len = max(len(str(int(y_min + i * (y_max - y_min) / num_ticks))) for i in range(num_ticks+1))
        x_label_offset = 40 + longest_x_len * 5
        y_label_offset = 30 + longest_y_len * 6

        # Legend sizing
        legend_font = 12
        swatch_size = 14
        swatch_gap = 10
        legend_row_gap = 6
        est_text_widths = [len(str(cat)) * 7 for cat in category_labels]
        legend_inner_width = max(swatch_size + swatch_gap + tw for tw in est_text_widths)

        # Ensure the legend title always fits — set a minimum width
        min_title_width = len("Legend") * 8 + 20  # rough estimate: 8px per char + padding
        legend_width = max(legend_inner_width + 20, min_title_width)

        legend_height = n_categories * (swatch_size + legend_row_gap) + 30

        # Plot area (centered, smaller footprint to leave room for labels/legend)
        plot_width = width * 0.65
        plot_height = height * 0.65
        plot_left = (width - plot_width) / 2
        plot_right = plot_left + plot_width
        plot_top = (height - plot_height) / 2
        plot_bottom = plot_top + plot_height

        # Group slot geometry with automatic bar sizing
        slot_width = (plot_right - plot_left) / n_groups
        group_gap_fraction = 0.2  # reserve 20% of slot for separation
        bar_area_width = slot_width * (1 - group_gap_fraction)
        bar_width_px = bar_area_width / n_categories
        total_bar_block = n_categories * bar_width_px

        # Legend position (outside plot, inside canvas)
        legend_x = plot_right + 40
        legend_y = plot_top

        buffer = io.StringIO()
        buffer.write(
            f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' "
            f"viewBox='0 0 {width} {height}' style='font-family:sans-serif'>\n"
        )
        buffer.write("  <rect width='100%' height='100%' fill='white'/>\n")

        # Accessibility
        buffer.write(f"<title>{title}</title>\n")
        buffer.write(f"<desc>Grouped bar chart with {n_groups} groups and {n_categories} categories</desc>\n")

        # Title
        buffer.write(
            f"<text x='{(plot_left+plot_right)/2}' y='{plot_top-20}' "
            f"text-anchor='middle' font-size='16'>{title}</text>\n"
        )

        # Axes
        buffer.write(f"<line x1='{plot_left}' y1='{plot_top}' x2='{plot_left}' y2='{plot_bottom}' stroke='black'/>\n")
        buffer.write(f"<line x1='{plot_left}' y1='{plot_bottom}' x2='{plot_right}' y2='{plot_bottom}' stroke='black'/>\n")

        # Y ticks
        for i in range(num_ticks + 1):
            val = y_min + i * (y_max - y_min) / num_ticks
            y = plot_bottom - (val - y_min) / (y_max - y_min) * (plot_bottom - plot_top)
            buffer.write(f"<line x1='{plot_left-5}' y1='{y}' x2='{plot_left}' y2='{y}' stroke='black'/>\n")
            buffer.write(f"<text x='{plot_left-10}' y='{y+4}' text-anchor='end' font-size='12'>{int(val)}</text>\n")

        # Axis labels
        buffer.write(
            f"<text x='{(plot_left+plot_right)/2}' y='{plot_bottom + x_label_offset}' "
            f"text-anchor='middle' font-size='14'>{x_label}</text>\n"
        )
        buffer.write(
            f"<text x='{plot_left - y_label_offset}' y='{(plot_top+plot_bottom)/2}' text-anchor='middle' "
            f"font-size='14' transform='rotate(-90 {plot_left - y_label_offset},{(plot_top+plot_bottom)/2})'>{y_label}</text>\n"
        )

        # Group labels
        for g, glabel in enumerate(group_labels):
            group_center = plot_left + (g + 0.5) * slot_width
            buffer.write(
                f"<text x='{group_center}' y='{plot_bottom + 20}' font-size='12' text-anchor='end' "
                f"transform='rotate(-45 {group_center},{plot_bottom + 20})'>{glabel}</text>\n"
            )

        # Bars
        for g, group_count in enumerate(counts):
            group_left = plot_left + g * slot_width
            start_x = group_left + (slot_width - bar_area_width) / 2
            group_total = sum(group_count.values()) or 1

            for c, cat in enumerate(categories):
                count = group_count.get(cat, 0)
                bar_h = (count - y_min) / (y_max - y_min) * (plot_bottom - plot_top)
                x_left = start_x + c * bar_width_px
                y_top = plot_bottom - bar_h

                buffer.write(
                    f"<rect x='{x_left}' y='{y_top}' width='{bar_width_px}' height='{bar_h}' "
                    f"fill='{colours[c]}' stroke='black'>\n"
                )
                frac = count / group_total
                buffer.write(f"  <title>{group_labels[g]} – {cat} – {count} ({frac:.1%})</title>\n")
                buffer.write("</rect>\n")
        # Legend
        buffer.write(
            f"<rect x='{legend_x}' y='{legend_y}' width='{legend_width}' height='{legend_height}' "
            f"fill='white' stroke='black' rx='6' ry='6'/>\n"
        )
        buffer.write(
            f"<text x='{legend_x + 10}' y='{legend_y + 18}' font-size='14' text-anchor='start'>Legend</text>\n"
        )
        for i, cat in enumerate(categories):
            # Add extra spacing after the legend title by including legend_row_gap
            entry_y = legend_y + 30 + legend_row_gap + i * (swatch_size + legend_row_gap)
            buffer.write(
                f"<rect x='{legend_x + 10}' y='{entry_y - swatch_size + 2}' "
                f"width='{swatch_size}' height='{swatch_size}' fill='{colours[i]}' stroke='black'/>\n"
            )
            buffer.write(
                f"<text x='{legend_x + 10 + swatch_size + swatch_gap}' y='{entry_y}' "
                f"font-size='{legend_font}' text-anchor='start'>{cat}</text>\n"
            )

        buffer.write("</svg>\n")

        # HTML wrapper with save button
        html_filename = filename
        html_content = f"""<!DOCTYPE html>
        <html>
        <head>
        <meta charset="utf-8">
        <title>{title}</title>
        <style>
            html, body {{
                height: 100%;
                margin: 0;
            }}
            .wrap {{
                width: 100vw;
                height: 100vh;
                display: flex;
                justify-content: center;
                align-items: center;
            }}
            svg {{
                max-width: 100vw;
                max-height: 100vh;
                width: 100%;
                height: auto;
                display: block;
                background: white;
                font-family: Arial, Helvetica, sans-serif;
            }}
            #saveBtn {{
                position: fixed;
                top: 10px;
                right: 10px;
                padding: 6px 12px;
                font-size: 14px;
                cursor: pointer;
            }}
        </style>
        </head>
        <body>
        <div class="wrap">
            {buffer.getvalue()}
        </div>
        <button id="saveBtn" onclick="downloadPNG()">Save as PNG</button>
        <script>
        function downloadPNG() {{
            const svg = document.querySelector("svg");
            const serializer = new XMLSerializer();
            const svgStr = serializer.serializeToString(svg);

            let width, height;
            if (svg.viewBox && svg.viewBox.baseVal) {{
                width = svg.viewBox.baseVal.width;
                height = svg.viewBox.baseVal.height;
            }} else {{
                width = svg.width.baseVal.value;
                height = svg.height.baseVal.value;
            }}

            const canvas = document.createElement("canvas");
            canvas.width = width;
            canvas.height = height;
            const ctx = canvas.getContext("2d");

            const img = new Image();
            img.onload = function() {{
                ctx.drawImage(img, 0, 0, width, height);
                const link = document.createElement("a");
                const safeTitle = "{title}".replace(/[^a-z0-9]+/gi, "_").toLowerCase();
                link.download = safeTitle + ".png";
                link.href = canvas.toDataURL("image/png");
                link.click();
            }};
            img.src = "data:image/svg+xml;base64," + btoa(unescape(encodeURIComponent(svgStr)));
        }}
        </script>
        </body>
        </html>"""

        with open(html_filename, "w", encoding="utf-8") as f:
            f.write(html_content)

        if open_browser:
            import webbrowser
            webbrowser.open(html_filename)

    def pie_chart(self,variable,title="Pie Chart",show_values=False,colours=None,legend=True,margin=60,filename="pie_chart.html",open_browser=False):
        """Create a pie chart from a list of categorical values ('variable') and configure chart elements.

        Parameters:

        variable (list): a list of categorical values.

        title (string): the title of the pie chart.

        show_values (Boolean): flag to indicate if counts and percentages should be shown in the legend beside each category label. If False, only category names are shown.

        colours (list): list of colour strings to use for slices. If None, colours are auto-generated. Colours may be specified as a common colour name or hex code.

        legend (Boolean): flag to indicate if a legend will be displayed.

        margin (integer): the margin size (in pixels) around the plot area.

        filename (string): the name of the output file. Include the output file path in the string or prepend tt.get_out_dir() to the filename.

        open_browser (Boolean): flag to indicate if the generated file should be automatically opened in the default web browser."""
        if not isinstance(variable, list) or not variable:
            raise ValueError("'variable' must be a non-empty list of categorical values.")
        if not isinstance(title, str):
            raise ValueError("'title' must be a string.")
        if not isinstance(show_values, bool):
            raise ValueError("'show_values' must be a Boolean value.")
        if colours is not None:
            if not isinstance(colours, list):
                raise ValueError("'colours' must be a list of colour strings.")
            if not all(isinstance(c, str) for c in colours):
                raise ValueError("All entries in 'colours' must be strings representing colours.")
        if not isinstance(legend, bool):
            raise ValueError("'legend' must be a Boolean value.")
        if not isinstance(margin, int) or margin < 0:
            raise ValueError("'margin' must be a non-negative integer.")
        if not isinstance(filename, str):
            raise ValueError("'filename' must be a string.")
        filename = self._normalize_html_filename(filename)
        if not isinstance(open_browser, bool):
            raise ValueError("'open_browser' must be a Boolean value.")
        from collections import Counter

        # Count categories
        counts = Counter(variable)
        categories = list(counts.keys())
        n_cat = len(categories)
        total = sum(counts.values())

        # Color palette handling
        def generate_hsl_palette(n):
            # evenly spaced hues with random offset, medium saturation/lightness for readable, distinct colours
            base_step = 360 / n
            offset = random.randint(0, 359)
            return [
                f"hsl({int((offset + i * base_step) % 360)}, 60%, 55%)"
                for i in range(n)
            ]

        if colours and len(colours) > 0:
            palette = colours[:]  # user-supplied
        else:
            palette = generate_hsl_palette(n_cat)

        # Dynamic margins (similar to line plot)
        extra_left = margin
        extra_right = margin + (100 if legend else 0)
        extra_top = margin + 40   # space for title
        extra_bottom = margin + 40  # space for legend if below

        width = 800 + extra_left + extra_right
        height = 600 + extra_top + extra_bottom

        plot_left = extra_left
        plot_right = width - extra_right
        plot_top = extra_top
        plot_bottom = height - extra_bottom

        # Pie geometry
        cx = (plot_left + plot_right) / 2
        cy = (plot_top + plot_bottom) / 2
        r = min(plot_right - plot_left, plot_bottom - plot_top) * 0.35

        # Begin SVG output
        buffer = io.StringIO()
        buffer.write(
            f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' "
            f"viewBox='0 0 {width} {height}' style='font-family:sans-serif'>\n"
        )
        buffer.write('  <rect width="100%" height="100%" fill="white"/>\n')

        # Title and description for accessibility
        buffer.write(f"<title>{title}</title>\n")
        buffer.write(f"<desc>Pie chart with {n_cat} categories</desc>\n")

        # Chart title
        buffer.write(
            f"<text x='{(plot_left+plot_right)/2}' y='{margin}' "
            f"text-anchor='middle' font-size='16'>{title}</text>\n"
        )

        # Draw slices
        start_angle = 0.0
        for i, category in enumerate(categories):
            count = counts[category]
            frac = count / total
            angle = frac * 2 * math.pi
            end_angle = start_angle + angle

            # Coordinates for arc
            x1 = cx + r * math.cos(start_angle)
            y1 = cy + r * math.sin(start_angle)
            x2 = cx + r * math.cos(end_angle)
            y2 = cy + r * math.sin(end_angle)
            large_arc = 1 if angle > math.pi else 0

            # Path for slice
            buffer.write(
                f"<path d='M {cx:.2f},{cy:.2f} L {x1:.2f},{y1:.2f} "
                f"A {r:.2f},{r:.2f} 0 {large_arc},1 {x2:.2f},{y2:.2f} Z' "
                f"fill='{palette[i % len(palette)]}'>\n"
            )
            # Tooltip
            buffer.write(
                f"  <title>{category} – {count} ({frac:.1%})</title>\n"
            )
            buffer.write("</path>\n")

            start_angle = end_angle

        # Legend (if enabled)
        if legend:
            legend_x = cx + r + 40   # closer to pie, inside canvas
            legend_y = plot_top + 20
            buffer.write("<g id='legend'>\n")
            for i, category in enumerate(categories):
                count = counts[category]
                frac = count / total
                label = str(category)
                if show_values:
                    label += f" – {count} ({frac:.1%})"

                y_offset = legend_y + i * 22
                buffer.write(
                    f"<rect x='{legend_x}' y='{y_offset}' width='16' height='16' "
                    f"fill='{palette[i % len(palette)]}'/>\n"
                )
                buffer.write(
                    f"<text x='{legend_x + 22}' y='{y_offset + 12}' font-size='14'>{label}</text>\n"
                )
            buffer.write("</g>\n")

        # Close SVG
        buffer.write("</svg>\n")

        # Also write HTML wrapper with save button
        html_filename = filename
        html_content = f"""<!DOCTYPE html>
        <html>
        <head>
        <meta charset="utf-8">
        <title>{title}</title>
        <style>
            html, body {{
                height: 100%;
                margin: 0;
            }}
            .wrap {{
                width: 100vw;
                height: 100vh;
                display: flex;
                justify-content: center;
                align-items: center;
            }}
            svg {{
                max-width: 100vw;
                max-height: 100vh;
                width: 100%;
                height: auto;
                display: block;
                background: white;
                font-family: Arial, Helvetica, sans-serif;
            }}
            #saveBtn {{
                position: fixed;
                top: 10px;
                right: 10px;
                padding: 6px 12px;
                font-size: 14px;
                cursor: pointer;
            }}
        </style>
        </head>
        <body>
        <div class="wrap">
            {buffer.getvalue()}
        </div>
        <button id="saveBtn" onclick="downloadPNG()">Save as PNG</button>
        <script>
        function downloadPNG() {{
            const svg = document.querySelector("svg");
            const serializer = new XMLSerializer();
            const svgStr = serializer.serializeToString(svg);

            let width, height;
            if (svg.viewBox && svg.viewBox.baseVal) {{
                width = svg.viewBox.baseVal.width;
                height = svg.viewBox.baseVal.height;
            }} else {{
                width = svg.width.baseVal.value;
                height = svg.height.baseVal.value;
            }}

            const canvas = document.createElement("canvas");
            canvas.width = width;
            canvas.height = height;
            const ctx = canvas.getContext("2d");

            const img = new Image();
            img.onload = function() {{
                ctx.drawImage(img, 0, 0, width, height);
                const link = document.createElement("a");
                const safeTitle = "{title}".replace(/[^a-z0-9]+/gi, "_").toLowerCase();
                link.download = safeTitle + ".png";
                link.href = canvas.toDataURL("image/png");
                link.click();
            }};
            img.src = "data:image/svg+xml;base64," + btoa(unescape(encodeURIComponent(svgStr)));
        }}
        </script>
        </body>
        </html>"""

        with open(html_filename, "w", encoding="utf-8") as f:
            f.write(html_content)

        # Optionally open in browser
        if open_browser:
            import webbrowser
            webbrowser.open(html_filename)

    # ============================================================
    # SPATIAL PLOTS
    # ============================================================

    def coordinate_plot(self, points, overlays = None, standard_dists = None, title="Coordinate Plot",x_label="X Coordinate",y_label="Y Coordinate",num_ticks=5,point_size=3,point_colours=None, overlay_colours = None, standard_dists_colours = None,legend=False,point_labels=None,overlay_labels = None,standard_dists_labels = None,tooltip_threshold=100,margin=60, filename="coordinate_plot.html",open_browser=False):
        """Create a coordinate plot for any number of PointDatasets ('points'), additional point overlays such as centroids ('overlays') and standard distances ("st) and configure plot elements.

        Parameters:

        points (object, list): the PointDataset or list of PointDatasets to be visualized. PointDatasets may contain any number of points. PointDatasets are the primary data layers to be visualized.

        overlays (object, tuple, list): the Point, list, tuple or list of point-like coordinate objects to be visualized.

        standard_dists (integer, float, list): the numerical value or list of numerical values representing the standard distances to be visualized. Standard distances (weighted and unweighted) can be calculated with the .points toolbox.

        title (string): the title of the coordinate plot. 

        x_label (string): the label of the x axis.

        y_label (string): the label of the y axis. 

        num_ticks (integer): the number of tick marks to draw on each axis.

        point_size (integer): the radius of points on the plot (in pixels).

        point_colours (string, list): if a single PointDataset is given, this should be a single string representing the colour of the points. If multiple PointDatasets are given, this should be a list of strings corresponding to the colour of each PointDataset. If no colours are specified or the number of specified colours does not equal the number of PointDatasets, random colours will be generated. Colours may be specified as a common colour name or hex code.
        
        overlay_colours (string, list): if a single point-like object is given, this should be a single string representing the colour of the points. If multiple point-like objects are given, this should be a list of strings corresponding to the colour of each point. If no colours are specified or the number of specified colours does not equal the number of points, random colours will be generated. Colours may be specified as a common colour name or hex code.
        
        standard_dists_colours (string, list): if a single standard distance is given, this should be a single string representing the colour of the standard distance to be visualized. If multiple standard distances are given, this should be a list of strings corresponding to the colour of each standard distance. If no colours are specified or the number of specified colours does not equal the number of standard distances, random colours will be generated. Colours may be specified as a common colour name or hex code.
        
        legend (Boolean): flag to indicate if a legend will be displayed.

        point_labels (string, list): if a single PointDataset is given, this should be a single string representing the label of the PointDataset. If multiple PointDatasets are given, this should be a list of strings representing the labels of each PointDataset. If no labels are given, generic labels will be generated. A legend must be created to view point labels.
        
        overlay_labels (string, list): if a single point-like object is given, this should be a single string representing the label of the point. If multiple point-like objects are given, this should be a list of strings representing the labels of each variable. If no labels are given, generic labels will be generated. A legend must be created to view point labels.
        
        standard_dists_labels (string, list): if a single standard distance is given, this should be a single string representing the label of the standard distance. If multiple standard distances are given, this should be a list of strings representing the labels of each standard distance. If no labels are given, generic labels will be generated. A legend must be created to view point labels.

        tooltip_threshold (integer): the maximum number of points per series for which explicit hover-based tooltips are included.
        
        margin (integer): the margin size (in pixels) around the plot area. 
        
        filename (string): the name of the output file. Include the output file path in the string or prepend tt.get_out_dir() to the filename.

        open_browser (Boolean): flag to indicate if the generated file should be automatically opened in the default web browser."""
        # Validate points
        if not isinstance(points, (list, tuple)):
            points = [points]
        if len(points) == 0:
            raise ValueError("'points' must contain at least one PointDataset.")
        for idx, pd in enumerate(points):
            if not isinstance(pd, _PointDataset):
                raise TypeError(f"All items in 'points' must be PointDataset objects.")
            if len(pd) == 0:
                raise ValueError(f"PointDataset at index {idx} is empty.")

        # Validate overlays
        if overlays is None:
            overlays = []
        else:
            overlays = [overlays]

        # overlays may contain: tuple, list, list of tuples, Point, PointDataset
        valid_overlay_types = (tuple, list, _Point, _PointDataset)
        for idx, ov in enumerate(overlays):
            if not isinstance(ov, valid_overlay_types):
                raise TypeError(
                    f"Overlay {idx} is not a valid point-like object. Overlays must be tuples, lists, lists of tuples, Point objects, or PointDatasets.")

        # Validate standard distances
        if standard_dists is None:
            standard_dists = []
        elif not isinstance(standard_dists, (list, tuple)):
            standard_dists = [standard_dists]

        for idx, sd in enumerate(standard_dists):
            if not isinstance(sd, (int, float)):
                raise TypeError(f"Standard distance {idx} must be a numeric value (int or float).")
            if sd <= 0:
                raise ValueError(f"Standard distance {idx} must be a positive numeric value.")

        # Validate title and axis labels
        if not isinstance(title, str):
            raise ValueError("'title' must be a string.")
        if not isinstance(x_label, str):
            raise ValueError("'x_label' must be a string.")
        if not isinstance(y_label, str):
            raise ValueError("'y_label' must be a string.")

        # Validate numeric parameters
        if not isinstance(num_ticks, int) or num_ticks <= 0:
            raise ValueError("'num_ticks' must be a positive integer.")
        if not isinstance(point_size, int) or point_size <= 0:
            raise ValueError("'point_size' must be a positive integer.")
        if not isinstance(tooltip_threshold, int) or tooltip_threshold < 0:
            raise ValueError("'tooltip_threshold' must be a non-negative integer.")
        if not isinstance(margin, int) or margin < 0:
            raise ValueError("'margin' must be a non-negative integer.")

        # Validate colours
        def _validate_colour_input(user_input, label):
            if user_input is None:
                return
            if isinstance(user_input, str):
                return
            if isinstance(user_input, (list, tuple)):
                if not all(isinstance(c, str) for c in user_input):
                    raise ValueError(f"All entries in '{label}' must be strings.")
                return
            raise ValueError(f"'{label}' must be a string or a list of strings.")

        _validate_colour_input(point_colours, "point_colours")
        _validate_colour_input(overlay_colours, "overlay_colours")
        _validate_colour_input(standard_dists_colours, "standard_dists_colours")

        # Validate legend
        if not isinstance(legend, bool):
            raise ValueError("'legend' must be a Boolean value.")

        # Validate labels
        def _validate_label_input(user_input, label):
            if user_input is None:
                return
            if isinstance(user_input, str):
                return
            if isinstance(user_input, (list, tuple)):
                if not all(isinstance(lbl, str) for lbl in user_input):
                    raise ValueError(f"All entries in '{label}' must be strings.")
                return
            raise ValueError(f"'{label}' must be a string or a list of strings.")

        _validate_label_input(point_labels, "point_labels")
        _validate_label_input(overlay_labels, "overlay_labels")
        _validate_label_input(standard_dists_labels, "standard_dists_labels")

        # Validate filename + browser flag
        if not isinstance(filename, str):
            raise ValueError("'filename' must be a string.")
        filename = self._normalize_html_filename(filename)
        if not isinstance(open_browser, bool):
            raise ValueError("'open_browser' must be a Boolean value.")

        # Collect all x and y coordinates from all PointDatasets
        all_x = []
        all_y = []

        for pd in points:
            for p in pd:
                all_x.append(p.x)
                all_y.append(p.y)

        # Collect overlay coordinates (normalize later, but we can extract raw coords now)
        def _extract_overlay_coords(obj):
            """Extract raw x,y pairs from any overlay object."""
            if isinstance(obj, _Point):
                return [(obj.x, obj.y)]
            elif isinstance(obj, _PointDataset):
                return [(p.x, p.y) for p in obj]
            elif isinstance(obj, tuple) and len(obj) == 2:
                return [obj]
            if isinstance(obj, list) and len(obj) == 2 and all(isinstance(v, (int, float)) for v in obj):
                return [(obj[0], obj[1])]
            elif isinstance(obj, list):
                # list of tuples or list of Points
                coords = []
                for item in obj:
                    if isinstance(item, _Point):
                        coords.append((item.x, item.y))
                    elif isinstance(item, tuple) and len(item) == 2:
                        coords.append(item)
                    elif isinstance(item, list) and len(item) == 2 and all(isinstance(v, (int, float)) for v in item):
                        coords.append((item[0], item[1]))
                    else:
                        raise TypeError(
                            "Overlay lists must contain only (x,y) tuples or Point objects."
                        )
                return coords
            else:
                raise TypeError("Invalid overlay object encountered during extent extraction.")

        for ov in overlays:
            coords = _extract_overlay_coords(ov)
            for x, y in coords:
                all_x.append(x)
                all_y.append(y)

        if len(all_x) == 0 or len(all_y) == 0:
            raise ValueError("No coordinate data found in 'points' or 'overlays'.")

        # Compute raw extents
        xmin = min(all_x)
        xmax = max(all_x)
        ymin = min(all_y)
        ymax = max(all_y)

        # Compute ranges
        x_range = xmax - xmin
        y_range = ymax - ymin

        # Avoid zero-range issues
        if x_range == 0:
            x_range = 1
        if y_range == 0:
            y_range = 1

        # Apply 10% buffer on each side
        buffer_x = x_range * 0.10
        buffer_y = y_range * 0.10

        xmin -= buffer_x
        xmax += buffer_x
        ymin -= buffer_y
        ymax += buffer_y

        # If standard distances are provided, expand extents to include circles
        # (We will normalize SDs later; for now assume SD is a radius and center is known)
        # NOTE: We only expand extents for datasets that have SDs.

        for idx, sd in enumerate(standard_dists):
            if idx >= len(points):
                break  # ignore extra SDs

            pd = points[idx]

            # Compute mean center for this dataset (same as your mean_center function)
            cx = sum(p.x for p in pd) / len(pd)
            cy = sum(p.y for p in pd) / len(pd)

            # Expand extents to include the SD circle
            xmin = min(xmin, cx - sd)
            xmax = max(xmax, cx + sd)
            ymin = min(ymin, cy - sd)
            ymax = max(ymax, cy + sd)
        
        # Base canvas size (consistent with other plot types)
        base_width, base_height = 800, 600

        # Estimate tick label width for dynamic Y-label offset
        tick_labels = [
            f"{ymin + i * (ymax - ymin) / num_ticks:.1f}"
            for i in range(num_ticks + 1)
        ]
        max_label_len = max(len(lbl) for lbl in tick_labels)
        font_size = 12  # matches tick label font size
        est_label_width = max_label_len * font_size * 0.6

        # Y-axis label offset (same pattern as scatter_plot)
        y_label_offset = margin - est_label_width - 15  # 15px padding

        # Extra left margin ensures rotated Y-label fits
        extra_left = max(0, margin - y_label_offset)

        # Extra right margin for legend (if enabled)
        extra_right = 180 if legend else 0

        # Final SVG dimensions
        width = base_width + margin * 2 + extra_left + extra_right
        height = base_height + margin * 2

        # Plot area boundaries
        plot_left   = margin + extra_left
        plot_right  = width - margin - extra_right
        plot_top    = margin
        plot_bottom = height - margin

        # Avoid division by zero (should already be handled, but double safety)
        if xmax == xmin:
            xmax = xmin + 1
        if ymax == ymin:
            ymax = ymin + 1

        # Linear scales
        xscale = (plot_right - plot_left) / (xmax - xmin)
        yscale = (plot_bottom - plot_top) / (ymax - ymin)

        # Coordinate transform functions
        sx = lambda x: plot_left + (x - xmin) * xscale
        sy = lambda y: plot_bottom - (y - ymin) * yscale

        num_datasets = len(points)
        num_overlays = len(overlays)
        num_sds = len(standard_dists)

        def random_colour_palette(n):
            """Generate n visually distinct HSL colours."""
            if n == 0:
                return []
            base_step = 360 / n
            offset = random.randint(0, 359)
            return [
                f"hsl({int((offset + i * base_step) % 360)},70%,50%)"
                for i in range(n)
            ]

        def resolve_colours(user_input, n, label):
            """
            TableTools colour rules:

            1. None → generate n random colours.
            2. String → only valid if n == 1.
            3. List/tuple → must contain exactly n strings.
            """

            # Case 1: No colours supplied
            if user_input is None:
                return random_colour_palette(n)

            # Case 2: Single string
            if isinstance(user_input, str):
                if n != 1:
                    raise ValueError(
                        f"'{label}' was given as a single string, but there are {n} items. "
                        f"Provide a list of {n} colours instead."
                    )
                return [user_input]

            # Case 3: List/tuple of strings
            if isinstance(user_input, (list, tuple)):
                if len(user_input) != n:
                    raise ValueError(
                        f"'{label}' must contain exactly {n} colours, "
                        f"but {len(user_input)} were provided."
                    )
                if not all(isinstance(c, str) for c in user_input):
                    raise ValueError(f"All entries in '{label}' must be strings.")
                return list(user_input)

            # Anything else is invalid
            raise ValueError(
                f"'{label}' must be None, a single string, or a list of {n} strings."
            )

        # Apply colour rules
        point_colours = resolve_colours(point_colours, num_datasets, "point_colours")
        overlay_colours = resolve_colours(overlay_colours, num_overlays, "overlay_colours")
        standard_dists_colours = resolve_colours(standard_dists_colours, num_sds, "standard_dists_colours")

        buf = io.StringIO()
        buf.write(
            f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}" '
            f'preserveAspectRatio="xMidYMid meet">\n'
        )
        buf.write('  <rect width="100%" height="100%" fill="white"/>\n')

        buf.write(f'  <title>{title}</title>\n')
        buf.write(f'  <desc>Coordinate plot with {len(points)} datasets, '
                f'{len(overlays)} overlays, and {len(standard_dists)} standard distances.</desc>\n')
        
        buf.write(
        f'  <text x="{width/2}" y="{margin/2}" '
        f'text-anchor="middle" font-size="20">{title}</text>\n')

        buf.write('  <defs>\n')

        # Each dataset gets its own marker shape (circle only for coordinate plot)
        # (We can extend this later if you want different shapes per dataset)
        for idx, colour in enumerate(point_colours):
            shape_def = f'<circle id="dataset_point{idx}" r="{point_size}" fill="{colour}"/>'
            buf.write(f'    {shape_def}\n')

        # Overlays get their own marker shape (slightly larger or different colour if desired)
        # For now: same size, but separate IDs so they can be styled independently
        for idx, colour in enumerate(overlay_colours):
            shape_def = f'<circle id="overlay_point{idx}" r="{point_size}" fill="{colour}"/>'
            buf.write(f'    {shape_def}\n')

        # Standard distance circles do not need <defs> entries (drawn directly)
        buf.write('  </defs>\n')

        # X-axis
        buf.write(
            f'  <line x1="{plot_left}" y1="{plot_bottom}" '
            f'x2="{plot_right}" y2="{plot_bottom}" stroke="black"/>\n'
        )

        # Y-axis
        buf.write(
            f'  <line x1="{plot_left}" y1="{plot_top}" '
            f'x2="{plot_left}" y2="{plot_bottom}" stroke="black"/>\n'
        )

        for i in range(num_ticks + 1):

            # X ticks
            tx = plot_left + i * (plot_right - plot_left) / num_ticks
            valx = xmin + i * (xmax - xmin) / num_ticks

            buf.write(
                f'  <line x1="{tx}" y1="{plot_bottom}" '
                f'x2="{tx}" y2="{plot_bottom+5}" stroke="black"/>\n'
            )
            buf.write(
                f'  <text x="{tx}" y="{plot_bottom+20}" '
                f'text-anchor="middle" font-size="12">{valx:.1f}</text>\n'
            )

            # Y ticks
            ty = plot_bottom - i * (plot_bottom - plot_top) / num_ticks
            valy = ymin + i * (ymax - ymin) / num_ticks

            buf.write(
                f'  <line x1="{plot_left}" y1="{ty}" '
                f'x2="{plot_left-5}" y2="{ty}" stroke="black"/>\n'
            )
            buf.write(
                f'  <text x="{plot_left-10}" y="{ty+4}" '
                f'text-anchor="end" font-size="12">{valy:.1f}</text>\n'
            )
        
        buf.write(
        f'  <text x="{(plot_left + plot_right)/2}" y="{plot_bottom + 45}" '
        f'text-anchor="middle" font-size="16">{x_label}</text>\n')

        y_label_x = plot_left - (extra_left - 10)
        y_label_y = (plot_top + plot_bottom) / 2

        buf.write(
            f'  <text x="{y_label_x}" y="{y_label_y}" '
            f'text-anchor="middle" font-size="16" '
            f'transform="rotate(-90,{y_label_x},{y_label_y})">{y_label}</text>\n'
        )

        if legend:

            legend_font = 12
            swatch_size = 14
            swatch_gap = 10
            legend_row_gap = 6

            # Build label lists (normalize later, but we need lengths now)
            num_datasets = len(points)
            num_overlays = len(overlays)
            num_sds = len(standard_dists)

            # Normalize labels (fallbacks)
            if point_labels is None:
                point_labels = [f"Dataset {i+1}" for i in range(num_datasets)]
            elif isinstance(point_labels,str):
                point_labels = [point_labels]
            else:
                point_labels = [
                    point_labels[i] if i < len(point_labels) else f"Dataset {i+1}"
                    for i in range(num_datasets)
                ]

            if overlay_labels is None:
                overlay_labels = [f"Overlay {i+1}" for i in range(num_overlays)]
            elif isinstance(overlay_labels,str):
                overlay_labels = [overlay_labels]
            else:
                overlay_labels = [
                    overlay_labels[i] if i < len(overlay_labels) else f"Overlay {i+1}"
                    for i in range(num_overlays)
                ]

            if standard_dists_labels is None:
                standard_dists_labels = [f"SD {i+1}" for i in range(num_sds)]
            elif isinstance(standard_dists_labels,str):
                standard_dists_labels = [standard_dists_labels]
            else:
                standard_dists_labels = [
                    standard_dists_labels[i] if i < len(standard_dists_labels) else f"SD {i+1}"
                    for i in range(num_sds)
                ]

            # Estimate text widths for all legend entries
            all_labels = point_labels + overlay_labels + standard_dists_labels
            est_text_widths = [len(lbl) * 7 for lbl in all_labels]

            legend_inner_width = max(swatch_size + swatch_gap + tw for tw in est_text_widths)

            # Ensure the "Legend" title always fits
            min_title_width = len("Legend") * 8 + 20
            legend_width = max(legend_inner_width + 20, min_title_width)

            # Total number of legend rows
            total_rows = num_datasets + num_overlays + num_sds
            legend_height = total_rows * (swatch_size + legend_row_gap) + 30

            # Legend position
            legend_x = plot_right + 40
            legend_y = plot_top

            # Legend frame
            buf.write(
                f"<rect x='{legend_x}' y='{legend_y}' width='{legend_width}' "
                f"height='{legend_height}' fill='white' stroke='black' rx='6' ry='6'/>\n"
            )
            buf.write(
                f"<text x='{legend_x + 10}' y='{legend_y + 18}' "
                f"font-size='14' text-anchor='start'>Legend</text>\n"
            )

            current_y = legend_y + 30 + legend_row_gap

            # Dataset entries
            for idx, label in enumerate(point_labels):
                colour = point_colours[idx]

                buf.write(
                    f"<circle cx='{legend_x + 10}' cy='{current_y}' "
                    f"r='{point_size}' fill='{colour}' stroke='black'/>\n"
                )
                buf.write(
                    f"<text x='{legend_x + 10 + swatch_size + swatch_gap}' "
                    f"y='{current_y+4}' font-size='{legend_font}' "
                    f"text-anchor='start'>{label}</text>\n"
                )

                current_y += swatch_size + legend_row_gap

            # Overlay entries
            for idx, label in enumerate(overlay_labels):
                colour = overlay_colours[idx]

                buf.write(
                    f"<circle cx='{legend_x + 10}' cy='{current_y}' "
                    f"r='{point_size}' fill='{colour}' stroke='black'/>\n"
                )
                buf.write(
                    f"<text x='{legend_x + 10 + swatch_size + swatch_gap}' "
                    f"y='{current_y+4}' font-size='{legend_font}' "
                    f"text-anchor='start'>{label}</text>\n"
                )

                current_y += swatch_size + legend_row_gap

            # Standard distance entries
            for idx, label in enumerate(standard_dists_labels):
                colour = standard_dists_colours[idx]

                # SD swatch: small circle outline
                buf.write(
                    f"<circle cx='{legend_x + 10}' cy='{current_y}' "
                    f"r='{point_size}' fill='none' stroke='{colour}' stroke-width='2'/>\n"
                )
                buf.write(
                    f"<text x='{legend_x + 10 + swatch_size + swatch_gap}' "
                    f"y='{current_y+4}' font-size='{legend_font}' "
                    f"text-anchor='start'>{label}</text>\n"
                )

                current_y += swatch_size + legend_row_gap

        for idx, pd in enumerate(points):
            buf.write(f'  <g id="dataset{idx}">\n')

            include_tooltips = len(pd) <= tooltip_threshold
            colour = point_colours[idx]

            if include_tooltips:
                # Explicit shapes with <title> tooltips
                for p in pd:
                    cx, cy = sx(p.x), sy(p.y)
                    buf.write(
                        f'    <circle cx="{cx}" cy="{cy}" r="{point_size}" '
                        f'fill="{colour}"><title>({p.x:.2f},{p.y:.2f})</title></circle>\n'
                    )
            else:
                # Compact <use> references (no tooltips)
                for p in pd:
                    buf.write(
                        f'    <use href="#dataset_point{idx}" '
                        f'x="{sx(p.x)}" y="{sy(p.y)}"/>\n'
                    )

            buf.write('  </g>\n')

        for idx, ov in enumerate(overlays):
            buf.write(f'  <g id="overlay{idx}">\n')

            colour = overlay_colours[idx]

            # Normalize overlay to a list of (x, y) pairs
            coords = _extract_overlay_coords(ov)

            include_tooltips = len(coords) <= tooltip_threshold

            if include_tooltips:
                for (x, y) in coords:
                    cx, cy = sx(x), sy(y)
                    buf.write(
                        f'    <circle cx="{cx}" cy="{cy}" r="{point_size}" '
                        f'fill="{colour}"><title>({x:.2f},{y:.2f})</title></circle>\n'
                    )
            else:
                for (x, y) in coords:
                    buf.write(
                        f'    <use href="#overlay_point{idx}" '
                        f'x="{sx(x)}" y="{sy(y)}"/>\n'
                    )

            buf.write('  </g>\n')

        for idx, sd in enumerate(standard_dists):

            if idx >= len(points):
                break  # ignore extra SDs

            pd = points[idx]
            colour = standard_dists_colours[idx]

            # Compute mean center (same logic as your mean_center function)
            cx = sum(p.x for p in pd) / len(pd)
            cy = sum(p.y for p in pd) / len(pd)

            # Convert center + radius to screen coordinates
            cx_screen = sx(cx)
            cy_screen = sy(cy)

            # Convert radius to screen units (x-scale only; circles remain circles)
            r_screen = sd * xscale

            buf.write(f'  <g id="standard_dist{idx}">\n')
            buf.write(
                f'    <circle cx="{cx_screen}" cy="{cy_screen}" r="{r_screen}" '
                f'fill="none" stroke="{colour}" stroke-width="2"/>\n'
            )
            buf.write('  </g>\n')    

        buf.write('</svg>\n')

        html_filename = filename
        html_content = f"""<!DOCTYPE html>
        <html>
        <head>
        <meta charset="utf-8">
        <title>{title}</title>
        <style>
            html, body {{
                height: 100%;
                margin: 0;
            }}
            .wrap {{
                width: 100vw;
                height: 100vh;
                display: flex;
                justify-content: center;
                align-items: center;
            }}
            svg {{
                max-width: 100vw;
                max-height: 100vh;
                width: 100%;
                height: auto;
                display: block;
                background: white;
                font-family: Arial, Helvetica, sans-serif;
            }}
            #saveBtn {{
                position: fixed;
                top: 10px;
                right: 10px;
                padding: 6px 12px;
                font-size: 14px;
                cursor: pointer;
            }}
        </style>
        </head>
        <body>
        <div class="wrap">
            {buf.getvalue()}
        </div>
        <button id="saveBtn" onclick="downloadPNG()">Save as PNG</button>
        <script>
            function downloadPNG() {{
                const svg = document.querySelector("svg");
                const serializer = new XMLSerializer();
                const svgStr = serializer.serializeToString(svg);

                // Use full viewBox dimensions if present
                let width, height;
                if (svg.viewBox && svg.viewBox.baseVal) {{
                    width = svg.viewBox.baseVal.width;
                    height = svg.viewBox.baseVal.height;
                }} else {{
                    width = svg.width.baseVal.value;
                    height = svg.height.baseVal.value;
                }}

                const canvas = document.createElement("canvas");
                canvas.width = width;
                canvas.height = height;
                const ctx = canvas.getContext("2d");

                const img = new Image();
                img.onload = function() {{
                    ctx.drawImage(img, 0, 0, width, height);
                    const link = document.createElement("a");
                    const safeTitle = "{title}".replace(/[^a-z0-9]+/gi, "_").toLowerCase();
                    link.download = safeTitle + ".png";
                    link.href = canvas.toDataURL("image/png");
                    link.click();
                }};
                img.src = "data:image/svg+xml;base64," + btoa(unescape(encodeURIComponent(svgStr)));
            }}
        </script>
        </body>
        </html>"""

        with open(html_filename, "w") as f:
            f.write(html_content)

        if open_browser:
            import webbrowser
            webbrowser.open(html_filename)

class _Matrix():
    """A set of functions for performing matrix operations and algebra on matrices."""
    # ============================================================
    # PRIVATE METHODS AND CLASS ATTRIBUTES AND FUNCTION DISPLAY
    # ============================================================
    def __init__(self):
        self._total_funcs = len([f for f in dir(__class__) if not f.startswith('__')])

    def __repr__(self):
        """If the toolbox is printed, display a message."""
        return ".matrix: a set of functions for performing matrix operations and algebra on matrices."

    def display_toolbox_functions(self):
        """Display a list of all available functions within this toolbox."""
        print(f"Number of {__class__.__name__[1:]} functions: {len([f for f in dir(__class__) if not f.startswith('__')])}")
        for f in [f for f in dir(__class__) if not f.startswith("__")]:
            print(f)
    
    # ============================================================
    # MATRIX CONVERSION
    # ============================================================
    
    def table_to_matrix(self,table):
        """Convert a Table object ('table') into a matrix and return a new Matrix object.
        
        Parameters:
        
        table (object): the Table object to be converted into a Matrix object."""
        if not isinstance(table, _Table):
            raise ValueError("Input must be a Table object.")

        return self.new_matrix(table.get_data())

    def matrix_to_table(self,matrix):
        """Convert a Matrix object ('matrix') into a Table object.
        
        Parameters:
        
        matrix (object): the Matrix object to be converted into a Table object."""
        if not isinstance(matrix, _MatrixObj):
            raise ValueError("Input must be a Matrix object.")
        output = TableTools().new_table()
        output.set_data(copy.deepcopy(matrix._matrix))
        output.autogenerate_headers()
        return output

    # ============================================================
    # MATRIX CREATION
    # ============================================================

    def new_matrix(self,rows):
        """Create a matrix from a nested list representing the rows of the matrix and return a Matrix object.
        
        Parameters:
        
        rows (list): a nested list representing the rows of the matrix."""
        if not isinstance(rows, list):
            raise ValueError("Input must be a list of rows.")
        if any(not isinstance(r, list) for r in rows):
            raise ValueError("Each row must be a list.")
        if len(set(len(r) for r in rows)) > 1:
            raise ValueError("All rows must have the same length.")

        m = _MatrixObj()
        m._matrix = copy.deepcopy(rows)
        m._initialize()
        return m
    
    def random_matrix(self,rows,columns,min,max):
        """Create a random Matrix object of the specified dimensions with values within the specified minimum ('min') and maximum ('max').
        
        Parameters:
        
        rows (integer): the number of rows in the output Matrix object.
        
        columns (integer): the number of columns in the output Matrix object.
        
        min (integer): the minimum value in the output Matrix object.
        
        max (integer): the maximum value in the output Matrix object."""
        if not isinstance(rows, int) or not isinstance(columns, int):
            raise ValueError("Rows and columns must be integers.")
        if rows <= 0 or columns <= 0:
            raise ValueError("Rows and columns must be positive.")
        if not isinstance(min, int) or not isinstance(max, int):
            raise ValueError("Min and max must be integers.")
        if min > max:
            raise ValueError("Minimum value cannot be greater than maximum value.")

        output = []
        if rows == 1:
            output = [random.randint(min, max) for _ in range(columns)]
        else:
            for _ in range(rows):
                output.append([random.randint(min, max) for _ in range(columns)])

        m = self.new_matrix(output)
        return m

    def identity(self,matrix):
        """Return the identity matrix of the input Matrix object ('matrix'). The identity matrix is a a matrix that, when multiplied by some matrix A, returns A.
        
        Parameters:
        
        matrix (object): the Matrix object the operation will be performed on. The input Matrix will not be modified by the operation."""
        if not isinstance(matrix, _MatrixObj):
            raise ValueError("Input must be a Matrix object.")
        if matrix.rows == 0 or matrix.columns == 0:
            raise ValueError("Cannot create identity matrix from empty matrix.")

        size = matrix.columns  # identity is based on number of columns
        output = []
        for i in range(size):
            row = [0 for _ in range(size)]
            row[i] = 1
            output.append(row)

        m = self.new_matrix(output)
        return m

    def zero_matrix(self,matrix):
        """Create a zero matrix to match the dimensions of an input Matrix object and return a new Matrix object.

        Parameters:
        
        matrix (object): the Matrix object the operation will be performed on. The input Matrix will not be modified by the operation."""
        if not isinstance(matrix, _MatrixObj):
            raise ValueError("Input must be a Matrix object.")
        if matrix.rows == 0 or matrix.columns == 0:
            raise ValueError("Cannot create zero matrix from empty dimensions.")

        output = [[0 for _ in range(matrix.columns)] for _ in range(matrix.rows)]
        return self.new_matrix(output)

    # ============================================================
    # ELEMENTARY ROW OPERATIONS
    # ============================================================

    def switch_rows(self,matrix,row_1,row_2):
        """Switch the position of any two rows in the input Matrix object ('matrix') and return a new Matrix object. Row switching is an elementary matrix row operation.
        
        Parameters:

        matrix (object): the Matrix object the operation will be performed on. The input Matrix will not be modified by the operation.
        
        row_1 (integer): the index of the first row to be switched.
        
        row_2 (integer): the index of the second row to be switched."""
        if not isinstance(matrix, _MatrixObj):
            raise ValueError("Input must be a Matrix object.")
        if not (0 <= row_1 < matrix.rows) or not (0 <= row_2 < matrix.rows):
            raise ValueError("Row indices are out of range.")

        # Deep copy to avoid mutating the original
        output = [row[:] for row in matrix._matrix]

        # Swap rows
        output[row_1], output[row_2] = output[row_2], output[row_1]

        return self.new_matrix(output)

    def row_multiplication(self,matrix,row,constant):
        """Multiply a specified row ('row') by a non-zero constant ('constant') and return a new Matrix object. Row multiplication is an elementary matrix row operation.
        
        Parameters:

        matrix (object): the Matrix object the operation will be performed on. The input Matrix will not be modified by the operation.
        
        row (integer): the index of the row to be multiplied.
        
        constant (integer, float): the non-zero constant to multiply the row values by."""
        if not isinstance(matrix, _MatrixObj):
            raise ValueError("Input must be a Matrix object.")
        if not (0 <= row < matrix.rows):
            raise ValueError("Row index is out of range.")
        if constant == 0:
            raise ValueError("Constant must not be zero.")

        # Deep copy to avoid mutating the original
        output = [r[:] for r in matrix._matrix]

        # Multiply the specified row
        output[row] = [val * constant for val in output[row]]

        return self.new_matrix(output)

    def row_addition(self,matrix,row_1,row_2):
        """Add the elements of two rows in a Matrix object ('matrix') and replace the first row. A new Matrix object is returned. Row addition is an elementary matrix row operation.
        
        Parameters:

        matrix (object): the Matrix object the operation will be performed on. The input Matrix will not be modified by the operation.
        
        row_1 (integer): the index of the first row to be added. This row will be replaced.
        
        row_2 (integer): the index of the second row to be added."""
        if not isinstance(matrix, _MatrixObj):
            raise ValueError("Input must be a Matrix object.")
        if not (0 <= row_1 < matrix.rows) or not (0 <= row_2 < matrix.rows):
            raise ValueError("Row indices are out of range.")

        r1 = matrix._as_row(row_1)
        r2 = matrix._as_row(row_2)
        if len(r1) != len(r2):
            raise ValueError("Rows must be of equal length.")

        o_row = [a + b for a, b in zip(r1, r2)]

        # Deep copy to avoid mutating the original
        output = [row[:] for row in matrix._matrix]
        output[row_1] = o_row

        return self.new_matrix(output)

    # ============================================================
    # SCALAR AND ELEMENTWISE OPERATIONS
    # ============================================================

    def matrix_addition(self,matrix_A,matrix_B):
        """Add the elements of two matrices element-wise and return a new Matrix object. If two matrices do not have equivalent dimensions, the addition is undefined and None will be returned.
        
        Parameters:

        matrix_A (object): the first Matrix object the operation will be performed on. The input Matrix will not be modified by the operation.

        matrix_B (object): the second Matrix object the operation will be performed on. The input Matrix will not be modified by the operation."""
        if not isinstance(matrix_A, _MatrixObj) or not isinstance(matrix_B, _MatrixObj):
            raise ValueError("Inputs must be Matrix objects.")
        if matrix_A.dimensions != matrix_B.dimensions:
            raise ValueError("Matrix dimensions must match for addition.")

        output = []
        for r in range(matrix_A.rows):
            s1 = matrix_A._as_row(r)
            o1 = matrix_B._as_row(r)
            o_row = [a + b for a, b in zip(s1, o1)]
            output.append(o_row)

        return self.new_matrix(output)
       
    def matrix_subtraction(self,matrix_A,matrix_B):
        """Subtract the elements of two matrices element-wise and return a new Matrix object. If two matrices do not have equivalent dimensions, the subtraction is undefined and None will be returned.
        
        Parameters:
        
        matrix_A (object): the first Matrix object the operation will be performed on. The input Matrix will not be modified by the operation.

        matrix_B (object): the second Matrix object the operation will be performed on. The input Matrix will not be modified by the operation."""
        if not isinstance(matrix_A, _MatrixObj) or not isinstance(matrix_B, _MatrixObj):
            raise ValueError("Inputs must be Matrix objects.")
        if matrix_A.dimensions != matrix_B.dimensions:
            raise ValueError("Matrix dimensions must match for subtraction.")

        output = []
        for r in range(matrix_A.rows):
            s1 = matrix_A._as_row(r)
            o1 = matrix_B._as_row(r)
            o_row = [a - b for a, b in zip(s1, o1)]
            output.append(o_row)

        return self.new_matrix(output)

    def scalar_multiplication(self,matrix,scalar):
        """Multiply each element of the input Matrix object by a scalar value and return a new Matrix object.
        
        Parameters:

        matrix (object): the Matrix object the operation will be performed on. The input Matrix will not be modified by the operation.
        
        scalar (integer, float): the scalar to multiply each element of the Matrix object by."""
        if not isinstance(matrix, _MatrixObj):
            raise ValueError("Input must be a Matrix object.")
        if not isinstance(scalar, (int, float)):
            raise ValueError("Scalar must be an integer or float.")

        output = [[matrix._matrix[row][col] * scalar for col in range(matrix.columns)]
                  for row in range(matrix.rows)]

        return self.new_matrix(output)

    # ============================================================
    # LINEAR ALGEBRA OPERATIONS
    # ============================================================

    def determinant(self,matrix):
        """Return the determinant of the input Matrix object ('matrix'). This function is currently only applicable to 2x2 matrices and 3x3 matrices, and None will be returned if a matrix with other dimensions is passed.
        
        Parameters:
        
        matrix (object): the Matrix object the operation will be performed on. The input Matrix will not be modified by the operation."""
        if not isinstance(matrix, _MatrixObj):
            raise ValueError("Input must be a Matrix object.")

        if matrix.dimensions == (2, 2):
            m = matrix._matrix
            a, b = m[0][0], m[0][1]
            c, d = m[1][0], m[1][1]
            return (a * d) - (b * c)

        elif matrix.dimensions == (3, 3):
            m = matrix._matrix
            a11, a12, a13 = m[0][0], m[0][1], m[0][2]
            a21, a22, a23 = m[1][0], m[1][1], m[1][2]
            a31, a32, a33 = m[2][0], m[2][1], m[2][2]
            return (a11 * a22 * a33) + (a21 * a32 * a13) + (a31 * a12 * a23) \
                   - (a11 * a32 * a23) - (a31 * a22 * a13) - (a21 * a12 * a33)

        else:
            raise ValueError("Determinant is only implemented for 2x2 and 3x3 matrices.")

    def inverse(self,matrix):
        """Return the inverse matrix of the input Matrix object ('matrix'). This function is currently only applicable to 2x2 matrices and 3x3 matrices, and None will be returned if a matrix with other dimensions is passed. If the determinant of the matrix is 0, the inverse is undefined and None will be returned.
        
        Parameters:
        
        matrix (object): the Matrix object the operation will be performed on. The input Matrix will not be modified by the operation."""
        if not isinstance(matrix, _MatrixObj):
            raise ValueError("Input must be a Matrix object.")

        if matrix.dimensions == (2, 2):
            m = matrix._matrix
            a, b = m[0][0], m[0][1]
            c, d = m[1][0], m[1][1]
            detM = (a * d) - (b * c)
            if detM == 0:
                return None
            adj = [[d, -b], [-c, a]]
            adj = self.new_matrix(adj)
            return self.scalar_multiplication(adj, 1 / detM)

        elif matrix.dimensions == (3, 3):
            m = matrix._matrix
            a11, a12, a13 = m[0][0], m[0][1], m[0][2]
            a21, a22, a23 = m[1][0], m[1][1], m[1][2]
            a31, a32, a33 = m[2][0], m[2][1], m[2][2]
            detM = (a11 * a22 * a33) + (a21 * a32 * a13) + (a31 * a12 * a23) \
                   - (a11 * a32 * a23) - (a31 * a22 * a13) - (a21 * a12 * a33)
            if detM == 0:
                return None
            adj = [
                [(a22 * a33) - (a23 * a32), (a13 * a32) - (a12 * a33), (a12 * a23) - (a13 * a22)],
                [(a23 * a31) - (a21 * a33), (a11 * a33) - (a13 * a31), (a13 * a21) - (a11 * a23)],
                [(a21 * a32) - (a22 * a31), (a12 * a31) - (a11 * a32), (a11 * a22) - (a12 * a21)]
            ]
            adj = self.new_matrix(adj)
            return self.scalar_multiplication(adj, 1 / detM)

        else:
            raise ValueError("Inverse is only implemented for 2x2 and 3x3 matrices.")

    def transpose(self,matrix):
        """Transpose the Matrix object ('matrix') and return a new Matrix object.
        
        Parameters:
        
        matrix (object): the Matrix object the operation will be performed on. The input Matrix will not be modified by the operation."""
        # transposition involves turning each row of the matrix in a column.
        if not isinstance(matrix, _MatrixObj):
            raise ValueError("Input must be a Matrix object.")
        if matrix.rows == 0 or matrix.columns == 0:
            raise ValueError("Cannot transpose an empty matrix.")

        if matrix.rows == 1:  # 1D row vector, column vector
            output = [[i] for i in matrix._as_row(0)]
        elif matrix.columns == 1:  # 1D column vector, row vector
            output = [matrix._as_col(0)]
        else:  # General case
            output = [[matrix._matrix[r][c] for r in range(matrix.rows)] for c in range(matrix.columns)]

        return self.new_matrix(output)

    def trace(self,matrix):
        """Calculate and return the summation of the diagonal elements of a square Matrix object ('matrix'). If the matrix object is not square (i.e., if the number of rows and columns are not equal), None will be returned.
        
        Parameters:

        matrix (object): the Matrix object the operation will be performed on. The input Matrix will not be modified by the operation."""
        if not isinstance(matrix, _MatrixObj):
            raise ValueError("Input must be a Matrix object.")
        if matrix.rows == 0 or matrix.columns == 0:
            raise ValueError("Cannot calculate trace of an empty matrix.")
        if matrix.rows != matrix.columns:
            raise ValueError("Trace is only defined for square matrices.")

        return sum(matrix._matrix[r][r] for r in range(matrix.rows))

    def inner_product(self,vector_A,vector_B):
        """Calculate the inner product of two vectors (1D matrices) and return a scalar value. If the vectors are not of equal length or vectors are not 1D, None will be returned.
        
        Parameters:
        
        vector_A (object): the first 1D matrix object the operation will be performed on. The input Matrix will not be modified by the operation.
    
        vector_B (object): the second 1D matrix object the operation will be performed on. The input Matrix will not be modified by the operation."""
        if not isinstance(vector_A, _MatrixObj) or not isinstance(vector_B, _MatrixObj):
            raise ValueError("Inputs must be Matrix objects.")

        # Ensure both are 1D
        if (vector_A.rows > 1 and vector_A.columns > 1) or (vector_B.rows > 1 and vector_B.columns > 1):
            raise ValueError("Inner product is only defined for 1D vectors.")

        # Convert column vectors to row vectors
        if vector_A.columns == 1:
            vector_A = self.transpose(vector_A)
        if vector_B.columns == 1:
            vector_B = self.transpose(vector_B)

        va = vector_A._matrix[0]
        vb = vector_B._matrix[0]

        if len(va) != len(vb):
            raise ValueError("Vectors must be of equal length.")

        return sum(map(operator.mul, va, vb))

    def matrix_multiplication(self, matrix_A, matrix_B):
        """Multiply two matrices ('matrix_A' and 'matrix_B') and return a new Matrix object. Matrix multiplication is only defined when the number of columns in 'matrix_A' equals the number of rows in 'matrix_B'. If the dimensions are incompatible, None will be returned.
        
        Parameters:
        
        matrix_A (object): the first Matrix object the operation will be performed on. The input Matrix will not be modified by the operation.
        
        matrix_B (object): the second Matrix object the operation will be performed on. The input Matrix will not be modified by the operation."""
        if not isinstance(matrix_A, _MatrixObj) or not isinstance(matrix_B, _MatrixObj):
            return None
        if matrix_A.columns != matrix_B.rows:
            return None

        output = []
        for i in range(matrix_A.rows):
            row = []
            for j in range(matrix_B.columns):
                val = sum(matrix_A._matrix[i][k] * matrix_B._matrix[k][j] for k in range(matrix_A.columns))
                row.append(val)
            output.append(row)

        return self.new_matrix(output)
      
# Matrix class for matrix operations and algebra
class _MatrixObj():
    """A matrix object to perform matrix operations and algebra with."""
    def __init__(self):
        self.rows = 0
        self.columns = 0
        self.dimensions = 0
        self._matrix = []
    
    def __repr__(self):
        """Print the matrix to the terminal when the matrix object is printed."""
        to_print = []
        out_str = "\n"
        for c in range(self.columns):
            col = self._as_col(c) #get each column as a list
            col = [str(v) for v in col] #convert each value to a string
            lvs = [len(v) for v in col] #find the length of each value, including the header
            ml = max(lvs) #find the max of the lengths of this column
            row = [col[i]+(" " * (ml-lvs[i]+4)) for i in range(len(col))] #add a number of blank spaces to the end of each value to equal the max column length + 2, so the columns are aligned
            to_print.append(row)
        for r in range(self.rows):
            row = []
            for c in range(self.columns):
                row.append(to_print[c][r]) #convert the corresponding values of each column into a row
            row = ''.join(row)
            out_str += row + "\n"
        return out_str
    
    def _initialize(self):
        """Initialize the matrix."""
        if isinstance(self._matrix[0],list): #if its 2D or a column vector
            self.rows = len(self._matrix)
            self.columns = len(self._matrix[0])
        else: # if it's a row vector
            self.rows = 1
            self.columns = len(self._matrix)
        self.dimensions = self.rows,self.columns
    
    def _as_row(self,row):
        """Get a row as a list."""
        if self.rows > 1 and self.columns > 1: #2D
            return self._matrix[row][:]
        elif self.rows > 1 and self.columns == 1: #1D column vector
            return self._matrix[row]
        elif self.rows == 1 and self.columns > 1: #1D row vector
            return self._matrix[:]
    
    def _as_col(self,col):
        """Get a column as a list."""
        output = []
        if self.rows > 1 and self.columns > 1: #2D
            for r in range(self.rows):
                output.append(self._matrix[r][col])
        elif self.rows > 1 and self.columns == 1: #1D column vector
            for r in range(self.rows):
                output.append(self._matrix[r][0])
        elif self.rows == 1 and self.columns > 1: #1D row vector
            output = [self._matrix[col]] 
        return output[:]

# Container class for point operations
class _PointDataset():
    """A container for a collection of Point objects, preserving table headers and dataset-level metadata (extent, schema mapping)."""
    def __init__(self, points, headers, id_col="id", x_col="x", y_col="y"):
        self.points = points
        self.headers = headers
        self.id_col = id_col
        self.x_col = x_col
        self.y_col = y_col
        xs = [p.x for p in points]
        ys = [p.y for p in points]
        self.extent = (min(xs), min(ys), max(xs), max(ys))
    
    def __repr__(self):
        """If the PointDataset object is returned, display a message."""
        return f"PointDataset Object: {len(self.points)} points."

    def __iter__(self):
        """Allow iteration directly over the points."""
        return iter(self.points)

    def __len__(self):
        """Allow returning the number of points directly."""
        return len(self.points)

    def __getitem__(self, idx):
        """Allow extraction of points directly."""
        return self.points[idx]
    
    def _update_extent(self):
        """Recalculate the dataset extent (min_x, min_y, max_x, max_y)."""
        if not self.points:
            self.extent = None
            return
        xs = [p.x for p in self.points]
        ys = [p.y for p in self.points]
        self.extent = (min(xs), min(ys), max(xs), max(ys))

# Individual point class for point operations
class _Point():
    """Represents a single point with enforced id, x, y attributes and a list of additional attribute values."""
    __slots__ = ("id", "x", "y", "attrs")
    def __init__(self,id,x,y,attrs=None):
        self.id = id
        self.x = x
        self.y = y
        self.attrs = attrs or []

    def __repr__(self):
        """If the Point object is returned, display a message."""
        return f"_Point(id={self.id}, x={self.x}, y={self.y}, attrs={self.attrs})"

    def __str__(self):
        """If the Point object is printed, display a message."""
        return self.__repr__()

# Symbol class for useful text symbols.
class _Symbols():
    """A symbol class object to hold useful text symbols to include in strings."""
    def __init__(self):
        self.superscript_zero = '⁰'
        self.superscript_one = '¹'
        self.superscript_two = '²'
        self.superscript_three = '³'
        self.superscript_four = '⁴'
        self.superscript_five = '⁵'
        self.superscript_six = '⁶'
        self.superscript_seven = '⁷'
        self.superscript_eight = '⁸'
        self.superscript_nine = '⁹'
        self.superscript_n = 'ⁿ'
        self.subscript_i = 'ᵢ'
        self.subscript_r = 'ᵣ'
        self.subscript_j = 'ⱼ'
        self.subscript_u = 'ᵤ'
        self.subscript_v = 'ᵥ'
        self.subscript_a = 'ₐ'
        self.subscript_e = 'ₑ'
        self.subscript_x = 'ₓ'
        self.degree = '°'
        self.plus_minus = '±'
        self.not_equal_to = '≠'
        self.less_than_equal_to = '≤'
        self.greater_than_equal_to = '≥'
        self.root = '√'
        self.infinity = '∞'
        self.alpha = 'α'
        self.beta = 'β'
        self.gamma = 'γ'
        self.delta = 'δ'
        self.delta_triangle = 'Δ'
        self.epsilon = 'ε'
        self.zeta = 'Ζ'
        self.eta = 'η'
        self.theta = 'θ'
        self.iota = 'ɩ'
        self.kappa = 'κ'
        self.lambda_sym = 'λ'
        self.mu = 'μ'
        self.nu = 'ν'
        self.xi = 'ξ'
        self.omicron = 'Ο'
        self.pi = 'π'
        self.rho = 'ρ'
        self.sigma = 'Σ'
        self.tau = '𝜏'
        self.upsilon = 'Ʊ'
        self.phi = 'φ'
        self.chi = 'χ'
        self.psi = 'ᴪ'
        self.omega = 'ω'
    
    def __repr__(self):
        """If the Symbols object is returned, display a message."""
        return f"TableTools symbol class"

# Table class
class _Table():
    """An object for storing and manipulating tabular data."""
    # ============================================================
    # PRIVATE METHODS AND CLASS ATTRIBUTES
    # ============================================================
    _SAFE_ENV = {
            "__builtins__": {},   # block unsafe builtins
            "abs": abs,
            "round": round,
            "min": min,
            "max": max,
            "len": len,
            "math": math,
            "str": str,
            "int": int,
            "float": float,
            "sum": sum,
            "bool": bool,
            "pow": pow,
            "any": any,
            "all": all
        } 
    # DICTIONARY FOR COLUMN CALCULATOR FUNCTION
    # - Function keys MUST be ALL CAPS.
    # - Values MUST be callables that accept a list of numeric values and return a single scalar.
    # - The calculator will call these functions exactly as written.
    _SCALAR_FUNCTIONS = {
        "HAS_DUPS": {
            "func": lambda col,: _ListOps().has_duplicates(col),
            "name": "Return True or False if the column has duplicates.",
            "args": 1
        },
        "SUM_SQR": {
            "func": lambda col, d: _Math().sum_squares(col,d),
            "name": "Sum of squared deviations from the mean.",
            "args": 1
        },
        "DIFF_SQR": {
            "func": lambda col1, col2, d: _Math().diff_squares(col1,col2,d),
            "name": "Difference of squared deviations between two columns.",
            "args": 2
        },
        "ACCUM": {
            "func": lambda col, d: _Math().accumulate(col,d),
            "name": "Cumulative sum",
            "args": 2
        },
        "MINAT": {
            "func": lambda col, thresh: _Math().min_above_threshold(col,thresh),
            "name": "Minimum value above a threshold.",
            "args": 2
        },
        "MAXAT": {
            "func": lambda col, thresh: _Math().max_above_threshold(col,thresh),
            "name": "Maximum value above a threshold.",
            "args": 2
        },
        "MINBT": {
            "func": lambda col, thresh: _Math().min_below_threshold(col,thresh),
            "name": "Minimum value below a threshold.",
            "args": 2
        },
        "MAXBT": {
            "func": lambda col, thresh: _Math().max_below_threshold(col,thresh),
            "name": "Maximum value below a threshold.",
            "args": 2
        },
        "FIRST": {
            "func": lambda col, d=None: col[0] if len(col) > 0 else None,
            "name": "First value in the column",
            "args": 1
        },

        "LAST": {
            "func": lambda col, d=None: col[-1] if len(col) > 0 else None,
            "name": "Last value in the column",
            "args": 1
        },
        "COUNT": {
            "func": lambda col, d: _Stats().descriptive_stat(col, "Count", d),
            "name": "Count of values",
            "args": 1
        },
        "UNIQUE": {
            "func": lambda col, d: _Stats().descriptive_stat(col, "Unique Values", d),
            "name": "Number of unique values",
            "args": 1
        },
        "SUM": {
            "func": lambda col, d: _Stats().descriptive_stat(col, "Sum", d),
            "name": "Sum of values",
            "args": 1
        },
        "MAX": {
            "func": lambda col, d: _Stats().descriptive_stat(col, "Max", d),
            "name": "Maximum value",
            "args": 1
        },
        "MIN": {
            "func": lambda col, d: _Stats().descriptive_stat(col, "Min", d),
            "name": "Minimum value",
            "args": 1
        },
        "RANGE": {
            "func": lambda col, d: _Stats().descriptive_stat(col, "Range", d),
            "name": "Range (max - min)",
            "args": 1
        },
        "MEAN": {
            "func": lambda col, d: _Stats().descriptive_stat(col, "Mean", d),
            "name": "Mean (average)",
            "args": 1
        },
        "MEDIAN": {
            "func": lambda col, d: _Stats().descriptive_stat(col, "Median", d),
            "name": "Median",
            "args": 1
        },
        "MODE": {
            "func": lambda col, d: _Stats().descriptive_stat(col, "Mode", d),
            "name": "Mode (most common value)",
            "args": 1
        },
        "IQ1": {
            "func": lambda col, d: _Stats().descriptive_stat(col, "IQ1", d),
            "name": "First quartile (Q1)",
            "args": 1
        },
        "IQ3": {
            "func": lambda col, d: _Stats().descriptive_stat(col, "IQ3", d),
            "name": "Third quartile (Q3)",
            "args": 1
        },
        "IQR": {
            "func": lambda col, d: _Stats().descriptive_stat(col, "IQR", d),
            "name": "Interquartile range (Q3 - Q1)",
            "args": 1
        },
        "STDEV": {
            "func": lambda col, d: _Stats().descriptive_stat(col, "Standard Deviation", d),
            "name": "Standard deviation",
            "args": 1
        },
        "STDERR": {
            "func": lambda col, d: _Stats().descriptive_stat(col, "Standard Error", d),
            "name": "Standard error of the mean",
            "args": 1
        },
        "VAR": {
            "func": lambda col, d: _Stats().descriptive_stat(col, "Variance", d),
            "name": "Variance",
            "args": 1
        },
        "COV": {
            "func": lambda col, d: _Stats().descriptive_stat(col, "Coefficient of Variation", d),
            "name": "Coefficient of variation",
            "args": 1
        },
        "SKEW": {
            "func": lambda col, d: _Stats().descriptive_stat(col, "Skewness", d),
            "name": "Skewness",
            "args": 1
        },
        "KURT": {
            "func": lambda col, d: _Stats().descriptive_stat(col, "Kurtosis", d),
            "name": "Kurtosis",
            "args": 1
        },
        "RMSE": {
            "func": lambda col1, col2, d: _Stats().root_mean_square_error(col1, col2, d),
            "name": "Root-mean-square-error",
            "args": 2
        },
        "MSE": {
            "func": lambda col1, col2, d: _Stats().mean_square_error(col1, col2, d),
            "name": "Mean-square-error",
            "args": 2
        },
        "MAE": {
            "func": lambda col1, col2, d: _Stats().mean_absolute_error(col1, col2, d),
            "name": "Mean-absolute-error",
            "args": 2
        },
        "FREQ": {
            "func": lambda col, value: _Stats().single_freq(col, value),
            "name": "Frequency of a single value.",
            "args": 2
        },
        "SPEARMAN": {
            "func": lambda col1, col2, d: _Stats().spearman_rank_correlation(col1, col2,d),
            "name": "Spearman Rank Correlation.",
            "args": 2
        },
    }

    # DICTIONARY FOR COLUMN CALCULATOR FUNCTION
    # - Column operation keys MUST be lower case.
    # - Values MUST be callables that accept a list of numeric values and return an equal-length list.
    # - The calculator will call these functions exactly as written.
    _COLUMN_FUNCTIONS = {
        "shuffle": {
            "func": lambda col: _ListOps().shuffle_list(col),
            "name": "Shuffle the order of values.",
            "args": 1
        },
        "cat_encode": {
            "func": lambda col, method: _ListOps().categorical_encoding(col, method),
            "name": "Encode categorical values.",
            "args": 2
        },
        "lin_rescale": {
            "func": lambda col, new_min, new_max, inverted, d: _Math().linear_rescale(col, new_min, new_max, inverted, d),
            "name": "Rescale values with a linear stretch.",
            "args": 4
        },
        "normalize": {
            "func": lambda col, d: _Math().normalize(col, d,),
            "name": "Normalize values from 0 to 1.",
            "args": 1
        },
        "min_max_clip": {
            "func": lambda col, min_val, max_val, d: _Math().min_max_clip(col, min_val, max_val, d,),
            "name": "Clip to the specified min and max values.",
            "args": 3
        },
        "percent_clip": {
            "func": lambda col, percent, tail, d: _Math().percent_clip(col, percent, tail, d,),
            "name": "Clip to the specified percentage and tail",
            "args": 3
        },
        "next_odd": {
            "func": lambda col, direction: _Math().next_odd(col, direction),
            "name": "Convert values to the next odd number in the specified direction",
            "args": 3
        },
        "next_even": {
            "func": lambda col, direction: _Math().next_even(col, direction),
            "name": "Convert values to the next even number in the specified direction",
            "args": 3
        },
        "even_odd": {
            "func": lambda col, mode, numeric: _Math().is_even_odd(col, mode, numeric),
            "name": "Evaluate if each numeric value is even or odd.",
            "args": 3
        },
        "cnvrt_dist": {
            "func": lambda col, from_unit, to_unit, d: _Conversion().convert_distance(col, from_unit, to_unit, d),
            "name": "Convert values between distance units.",
            "args": 3
        },
        "cnvrt_area": {
            "func": lambda col, from_unit, to_unit, d: _Conversion().convert_area(col, from_unit, to_unit, d),
            "name": "Convert values between area units.",
            "args": 3
        },
        "cnvrt_vol": {
            "func": lambda col, from_unit, to_unit, d: _Conversion().convert_volume(col, from_unit, to_unit, d),
            "name": "Convert values between volume units.",
            "args": 3
        },
        "cnvrt_wgt": {
            "func": lambda col, from_unit, to_unit, d: _Conversion().convert_weight(col, from_unit, to_unit, d),
            "name": "Convert values between weight units.",
            "args": 3
        },
        "cnvrt_dnst": {
            "func": lambda col, from_unit, to_unit, d: _Conversion().convert_density(col, from_unit, to_unit, d),
            "name": "Convert values between density units.",
            "args": 3
        },
        "cnvrt_time": {
            "func": lambda col, from_unit, to_unit, d: _Conversion().convert_time(col, from_unit, to_unit, d),
            "name": "Convert values between time units.",
            "args": 3
        },
        "cnvrt_vel": {
            "func": lambda col, from_unit, to_unit, d: _Conversion().convert_velocity(col, from_unit, to_unit, d),
            "name": "Convert values between velocity units.",
            "args": 3
        },
        "cnvrt_flow": {
            "func": lambda col, from_unit, to_unit, d: _Conversion().convert_flowrate(col, from_unit, to_unit, d),
            "name": "Convert values between flowrate units.",
            "args": 3
        },
        "cnvrt_prsr": {
            "func": lambda col, from_unit, to_unit, d: _Conversion().convert_pressure(col, from_unit, to_unit, d),
            "name": "Convert values between pressure units.",
            "args": 3
        },
        "cnvrt_enrg": {
            "func": lambda col, from_unit, to_unit, d: _Conversion().convert_energy(col, from_unit, to_unit, d),
            "name": "Convert values between energy units.",
            "args": 3
        },
        "cnvrt_pow": {
            "func": lambda col, from_unit, to_unit, d: _Conversion().convert_power(col, from_unit, to_unit, d),
            "name": "Convert values between power units.",
            "args": 3
        },
        "cnvrt_frc": {
            "func": lambda col, from_unit, to_unit, d: _Conversion().convert_force(col, from_unit, to_unit, d),
            "name": "Convert values between force units.",
            "args": 3
        },
        "cnvrt_temp": {
            "func": lambda col, from_unit, to_unit, d: _Conversion().convert_temperature(col, from_unit, to_unit, d),
            "name": "Convert values between temperature units.",
            "args": 3
        },
        "cnvrt_conc": {
            "func": lambda col, from_unit, to_unit, d: _Conversion().convert_concentration(col, from_unit, to_unit, d),
            "name": "Convert values between concentration units.",
            "args": 3
        },
        "standardize": {
            "func": lambda col, d: _Stats().standardize(col, d),
            "name": "Standard score/z-score.",
            "args": 1
        },
        "case": {
            "func": lambda col, mode: _Text().case(col, mode),
            "name": "Convert text case.",
            "args": 2
        },
        "capitalize": {
            "func": lambda col: _Text().capitalize(col,),
            "name": "Capitalize text.",
            "args": 1
        },
        "reverse": {
            "func": lambda col: _Text().reverse(col,),
            "name": "Reverse text.",
            "args": 1
        },
        "trunc": {
            "func": lambda col, length, ell: _Text().truncate(col,length, ell),
            "name": "Truncate text.",
            "args": 3
        },
        "char_dex": {
            "func": lambda col, dex: _Text().char_at_index(col,dex),
            "name": "The character at the specified index.",
            "args": 2
        },
        "dex_char": {
            "func": lambda col, char: _Text().index_of_char(col,char),
            "name": "The index of the specified character.",
            "args": 2
        },
        "substr": {
            "func": lambda col, start, end: _Text().extract_substring(col,start, end),
            "name": "Extract substring at specified indices.",
            "args": 3
        },
        "str_pos": {
            "func": lambda col, substr, mode: _Text().match_position(col,substr, mode),
            "name": "Evaluate if texts startswith/endwith specified substring.",
            "args": 3
        },
        "length": {
            "func": lambda col : _Text().length(col),
            "name": "Length of text.",
            "args": 1
        },
        "count_str": {
            "func": lambda col, substr : _Text().count_str(col, substr),
            "name": "Frequency of a substring.",
            "args": 2
        },
        "contains": {
            "func": lambda col, substr : _Text().contains(col, substr),
            "name": "Evaluate if text contains specified substring.",
            "args": 2
        },
        "replace": {
            "func": lambda col, old_str, new_str : _Text().replace(col, old_str, new_str),
            "name": "Replace a substring with another substring.",
            "args": 3
        },
        "whitespace": {
            "func": lambda col : _Text().remove_whitespace(col),
            "name": "Remove all whitespace characters from text.",
            "args": 1
        },
        "pad": {
            "func": lambda col, length, char, mode : _Text().pad(col, length, char, mode),
            "name": "Pad text to the specified length with a specified character.",
            "args": 4
        },
        "trim": {
            "func": lambda col, chars, mode : _Text().trim(col, chars, mode),
            "name": "Trim the specified number of characters from text.",
            "args": 3
        },
        "norm_date": {
            "func": lambda col, sep : _Date().normalize_dates(col, sep),
            "name": "Normalize date strings to ISO.",
            "args": 3
        },
    }

    def __init__(self, data = None):
        self._headers = []
        self._dtypes = []
        self._data = self._initialize_data(data)
        self._boolean_true_values = None
        self._boolean_false_values = None
        self._total_funcs = len([f for f in dir(__class__) if not f.startswith('_')])

    def __repr__(self):
        """If the Table object is returned, display a message."""
        return "TableTools Table Object"

    def __getitem__(self, key):
        """Allow users access to data via keys and route to the correct function."""
        if isinstance(key, int):
            return self.get_row(key)

        if isinstance(key, str):
            return self.get_column(key)

        if isinstance(key, tuple) and len(key) == 2:
            row, col = key
            return self.get_value(row, col)

        if isinstance(key, slice):
            return [self.get_row(i) for i in range(*key.indices(self.num_rows()))]

        raise TypeError("Invalid index type for Table.")

    def __setitem__(self, key, value):
        """Direct assignment via indexing is not supported. This restriction prevents accidental structural changes, dtype corruption, and bypassing of TableTools' explicit mutation methods."""
        raise TypeError("Direct assignment via indexing is not supported. Use existing methods for appending, inserting, and replacing rows/columns/values.")

    def _load_on_read(self,data):
        """Load data to a new Table object when reading data from a file.
        
        Parameters:
        
        data (list): a nested list of rows."""
        self._data = copy.deepcopy(data)

    def _initialize_data(self,data):
        """Initialize the data in a new Table object when not read from a file.
        
        Parameters:
        
        data (object): the Table object to extract data from."""
        if data:
            if hasattr(data, "get_headers") and hasattr(data, "get_dtypes") and hasattr(data, "get_data"):
                self._headers = data.get_headers()
                self._dtypes = data.get_dtypes()
                return data.get_data()[:]
            else:
                raise TypeError("Expected a Table object or None")
        return []
  
    def _return_col_type(self,col):
        """Return the data type of the column. In a column of mixed numerical values (floats and integers), all numerical values will be converted to floats.
        
        Parameters:
        
        col (list): a list of values to determine the data type of."""
        if not col:  # empty list or None
            raise ValueError("Column is empty. Cannot determine data type.")
        if any([isinstance(v,float) for v in col]):
            return "float"
        elif any(isinstance(v, int) and not isinstance(v, bool) for v in col):
            return "integer"
        elif any(isinstance(v, bool) for v in col):
            return "bool"
        elif any([isinstance(v,str) for v in col]):
            return "string"
  
    def _row_dtypes(self,vals):
        """Determine the data type of each value in a row and set the Table data types. Used when adding a first row of data to an empty Table.
        
        Parameters:
        
        vals (list): the list of values to determine data types for."""
        for v in vals:
            if isinstance(v, bool):
                self._dtypes.append("bool")
            elif isinstance(v, int):
                self._dtypes.append("integer")
            elif isinstance(v, float):
                self._dtypes.append("float")
            else:
                # everything else becomes string
                self._dtypes.append("string")

    def _viewer_html(self):
        """Produce the HTML code that the interactive viewer depends on."""
        return """<!doctype html>
        <html lang="en">
        <head>
        <meta charset="utf-8">
        <title>TableTools Viewer</title>
        <meta name="viewport" content="width=device-width, initial-scale=1">

        <style>
            body {
                font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
                margin: 0;
                background-color: #f4f7f6;
                color: #333;
            }
            header {
                padding: 12px 16px;
                display: flex;
                align-items: center;
                gap: 10px;
                position: sticky;
                top: 0;
                z-index: 10;
                background-color: #2c3e50;
                border-bottom: 2px solid #3498db;
                color: #ecf0f1;

            }
            header strong {
                font-size: 16px;
                color: #ecf0f1;
            }
            button {
                padding: 6px 12px;
                border-radius: 6px;
                cursor: pointer;
                font-size: 14px;
                background-color: #3498db;
                border: 1px solid #2980b9;
                color: #ecf0f1;

            }
            button:hover {
                background-color: #2980b9;
            }
            button:focus {
                outline: 2px solid #ecf0f1;
                outline-offset: 1px;
            }
            #grid {
                height: calc(100vh - 60px);
                overflow: auto;
                background: #ffffff;
                position: relative;
                white-space: nowrap;
            }
            table {
                border-collapse: collapse;
                font-size: 14px;
                table-layout: fixed;
                width: max-content; 
                min-width: 100%
            }
            th, td {
                border-bottom: 1px solid #eee;
                padding: 6px 10px;
                white-space: nowrap;
                overflow: hidden;
                text-overflow: ellipsis;
            }
            thead th {
                position: sticky;
                top: 0;
                background: #f0f0f0;
                border-bottom: 2px solid #ccc;
                z-index: 5;
                text-align: left;
                position: relative;
            }
            tbody tr:hover td {
                background: #eef6ff;
            }
            .selected-cell {
                outline: 2px solid #4a90e2;
                outline: 2px solid #3498db;
                background-color: #ecf0f1
            }
            .editing-cell {
                outline-offset: -2px;
                background-color: #ffffff !important;
                outline: 2px solid #2980b9
            }
            #spacer {
                position: absolute;
                left: 0;
                width: 1px;
                background: transparent;
                z-index: -1;
            }
            .col-resize-handle {
                position: absolute;
                top: 0;
                right: 0;
                width: 6px;
                height: 100%;
                cursor: col-resize;
                z-index: 10;
            }
            .col-resize-handle:hover {
                background: rgba(0, 0, 0, 0.05);
                background-color: rgba(52, 73, 94, 0.2)
            }
        </style>
        </head>

        <body>
        <header>
            <strong>TableTools Viewer</strong>
            <button id="saveBtn" aria-label="Save changes to the table">Save</button>
        </header>

        <div id="grid" role="region" aria-label="Editable data table">
            <div id="spacer"></div>
            <table id="table">
                <colgroup id="colgroup"></colgroup>
                <thead id="thead"></thead>
                <tbody id="tbody"></tbody>
            </table>
        </div>

        <script>
        let tableData = null;
        let edits = [];
        const rowHeight = 28;
        const buffer = 20;
        let totalRows = 0;
        let totalCols = 0;

        // Selection / editing state
        let selectedRow = 0;
        let selectedCol = 0;
        let editing = false;
        let editingCell = null;
        let originalEditValue = '';

        function escapeHtml(s) {
            return s.replace(/[&<>"']/g, m => ({
                '&':'&amp;', '<':'&lt;', '>':'&gt;', '"':'&quot;', "'":'&#39;'
            }[m]));
        }

        async function loadData() {
            const res = await fetch('/data');
            tableData = await res.json();
            totalRows = tableData.rows.length;
            totalCols = tableData.headers.length;

            // Clamp selection if table is empty
            if (totalRows === 0 || totalCols === 0) {
                selectedRow = 0;
                selectedCol = 0;
            } else {
                selectedRow = 0;
                selectedCol = 0;
            }

            renderTableShell();
            updateVisibleRows();
        }

        function renderTableShell() {
            const headers = tableData.headers;
            const thead = document.getElementById('thead');
            const colgroup = document.getElementById('colgroup');

            // Build colgroup for column resizing
            let colsHtml = '';
            for (let i = 0; i < headers.length; i++) {
                colsHtml += '<col>';
            }
            colgroup.innerHTML = colsHtml;

            // Build header row with resize handles
            let h = '<tr>';
            for (let i = 0; i < headers.length; i++) {
                const col = headers[i];
                h += (
                    '<th scope="col" data-col-index="' + i + '">' +
                    escapeHtml(String(col)) +
                    '<div class="col-resize-handle" data-col-index="' + i + '"></div>' +
                    '</th>'
                );
            }
            h += '</tr>';
            thead.innerHTML = h;

            const spacer = document.getElementById('spacer');
            spacer.style.height = (totalRows * rowHeight) + 'px';

            attachResizeHandlers();
        }

        function updateVisibleRows() {
            const grid = document.getElementById('grid');
            const tbody = document.getElementById('tbody');

            const scrollTop = grid.scrollTop;
            const viewportHeight = grid.clientHeight;

            const start = Math.max(0, Math.floor(scrollTop / rowHeight) - buffer);
            const end = Math.min(
                totalRows,
                Math.floor((scrollTop + viewportHeight) / rowHeight) + buffer
            );

            let html = '';
            for (let r = start; r < end; r++) {
                html += '<tr>';
                const row = tableData.rows[r];
                for (let c = 0; c < row.length; c++) {
                    const val = row[c] == null ? '' : String(row[c]);

                    let classes = '';
                    let contentEditable = 'false';
                    let tabIndex = '-1';

                    if (!editing && r === selectedRow && c === selectedCol) {
                        classes += ' selected-cell';
                        tabIndex = '0';
                    }

                    if (editing && r === selectedRow && c === selectedCol) {
                        classes += ' editing-cell';
                        contentEditable = 'true';
                        tabIndex = '0';
                    }

                    html += (
                        '<td ' +
                        'data-row="' + r + '" data-col="' + c + '" ' +
                        'class="' + classes.trim() + '" ' +
                        'contenteditable="' + contentEditable + '" ' +
                        'tabindex="' + tabIndex + '">' +
                        escapeHtml(val) + '</td>'
                    );
                }
                html += '</tr>';
            }

            tbody.innerHTML = html;
            tbody.style.transform = 'translateY(' + (start * rowHeight) + 'px)';

            attachCellHandlers();
        }

        function attachCellHandlers() {
            const cells = document.querySelectorAll('#tbody td');

            cells.forEach(cell => {
                const r = parseInt(cell.dataset.row, 10);
                const c = parseInt(cell.dataset.col, 10);

                cell.addEventListener('click', () => {
                    if (editing) {
                        // Finish current edit before changing selection
                        finishEdit(true);
                    }
                    setSelection(r, c);
                });

                cell.addEventListener('dblclick', () => {
                    if (!editing) {
                        setSelection(r, c);
                        startEdit(false, null);
                    }
                });
            });
        }

        function setSelection(r, c) {
            if (r < 0 || r >= totalRows || c < 0 || c >= totalCols) return;
            selectedRow = r;
            selectedCol = c;
            editing = false;
            editingCell = null;
            updateVisibleRows();
        }

        function focusSelectedCellIfVisible() {
            const selector = '#tbody td[data-row="' + selectedRow + '"][data-col="' + selectedCol + '"]';
            const cell = document.querySelector(selector);
            if (cell) {
                cell.focus();
            }
        }

        function ensureCellVisible() {
            const grid = document.getElementById('grid');
            if (!grid) return;

            // Compute the vertical position of the selected row
            const targetTop = selectedRow * rowHeight;
            const targetBottom = targetTop + rowHeight;

            const scrollTop = grid.scrollTop;
            const viewportHeight = grid.clientHeight;

            // Vertical scroll adjustment
            if (targetTop < scrollTop) {
                grid.scrollTop = targetTop;
            } else if (targetBottom > scrollTop + viewportHeight) {
                grid.scrollTop = targetBottom - viewportHeight;
            }

            // Horizontal scroll adjustment (based on column widths)
            const colgroup = document.getElementById('colgroup');
            let targetLeft = 0;
            for (let i = 0; i < selectedCol; i++) {
                const col = colgroup.children[i];
                const w = col && col.style.width ? parseInt(col.style.width) : 100;
                targetLeft += w;
            }
            const targetRight = targetLeft + 100;

            const scrollLeft = grid.scrollLeft;
            const viewportWidth = grid.clientWidth;

            if (targetLeft < scrollLeft) {
                grid.scrollLeft = targetLeft;
            } else if (targetRight > scrollLeft + viewportWidth) {
                grid.scrollLeft = targetRight - viewportWidth;
            }

            // Now that we've scrolled, re-render rows so the cell exists
            updateVisibleRows();

            // Now the cell exists — focus it
            const selector = '#tbody td[data-row="' + selectedRow + '"][data-col="' + selectedCol + '"]';
            const cell = document.querySelector(selector);
            if (cell) cell.focus();
        }

        // Editing logic

        function startEdit(replaceContent, initialChar) {
            if (totalRows === 0 || totalCols === 0) return;

            const selector = '#tbody td[data-row="' + selectedRow + '"][data-col="' + selectedCol + '"]';
            const cell = document.querySelector(selector);
            if (!cell) return; // Not in view yet

            editing = true;
            editingCell = cell;
            originalEditValue = cell.textContent;

            cell.classList.add('editing-cell');
            cell.classList.remove('selected-cell');
            cell.setAttribute('contenteditable', 'true');

            // If we want to replace content when typing starts
            if (replaceContent) {
                cell.textContent = initialChar != null ? initialChar : '';
            }

            cell.focus();

            // Move caret to end
            const range = document.createRange();
            range.selectNodeContents(cell);
            range.collapse(false);
            const sel = window.getSelection();
            sel.removeAllRanges();
            sel.addRange(range);

            // Attach key handling specific to edit mode
            cell.addEventListener('keydown', editModeKeyHandler);
        }

        function editModeKeyHandler(e) {
            if (!editing || !editingCell) return;

            const isShift = e.shiftKey;

            if (e.key === 'Enter') {
                e.preventDefault();
                finishEdit(true);
                // Move selection down or up depending on Shift
                const newRow = isShift ? selectedRow - 1 : selectedRow + 1;
                if (newRow >= 0 && newRow < totalRows) {
                    setSelection(newRow, selectedCol);
                } else {
                    setSelection(selectedRow, selectedCol);
                }
                return;
            }

            if (e.key === 'Tab') {
                e.preventDefault();
                finishEdit(true);
                const delta = isShift ? -1 : 1;
                let newRow = selectedRow;
                let newCol = selectedCol + delta;

                if (newCol < 0) {
                    if (newRow > 0) {
                        newRow -= 1;
                        newCol = totalCols - 1;
                    } else {
                        newCol = 0;
                    }
                } else if (newCol >= totalCols) {
                    if (newRow < totalRows - 1) {
                        newRow += 1;
                        newCol = 0;
                    } else {
                        newCol = totalCols - 1;
                    }
                }
                setSelection(newRow, newCol);
                return;
            }

            if (e.key === 'Escape') {
                e.preventDefault();
                // Cancel changes and restore original
                if (editingCell) {
                    editingCell.textContent = originalEditValue;
                }
                finishEdit(false);
                setSelection(selectedRow, selectedCol);
                return;
            }

            // Arrow keys in edit mode are left to the browser to move caret
        }

        function finishEdit(commit) {
            if (!editing || !editingCell) return;

            editingCell.removeEventListener('keydown', editModeKeyHandler);

            const r = selectedRow;
            const c = selectedCol;
            const newVal = editingCell.textContent;
            const oldVal = tableData.rows[r][c] == null ? '' : String(tableData.rows[r][c]);

            if (commit && newVal !== oldVal) {
                tableData.rows[r][c] = newVal;
                // Replace any existing edit for this cell
                edits = edits.filter(e => !(e.row === r && e.col === c));
                edits.push({ row: r, col: c, value: newVal });
            } else {
                // Revert visual content to model if not committing
                editingCell.textContent = oldVal;
            }

            editingCell.classList.remove('editing-cell');
            editingCell.setAttribute('contenteditable', 'false');
            editingCell = null;
            editing = false;

            updateVisibleRows();
        }

        // Navigation in selection mode

        document.addEventListener('keydown', (e) => {
            if (editing) {
                // Edit mode keys are handled in editModeKeyHandler
                return;
            }

            if (totalRows === 0 || totalCols === 0) return;

            const key = e.key;
            const isShift = e.shiftKey;
            const isMeta = e.metaKey || e.ctrlKey || e.altKey;

            // Typing in selection mode: start editing and replace content
            if (!isMeta && key.length === 1) {
                e.preventDefault();
                startEdit(true, key);
                return;
            }

            if (key === 'F2') {
                e.preventDefault();
                startEdit(false, null);
                return;
            }

            let handled = false;
            let newRow = selectedRow;
            let newCol = selectedCol;

            if (key === 'ArrowUp') {
                newRow = Math.max(0, selectedRow - 1);
                handled = true;
            } else if (key === 'ArrowDown') {
                newRow = Math.min(totalRows - 1, selectedRow + 1);
                handled = true;
            } else if (key === 'ArrowLeft') {
                newCol = Math.max(0, selectedCol - 1);
                handled = true;
            } else if (key === 'ArrowRight') {
                newCol = Math.min(totalCols - 1, selectedCol + 1);
                handled = true;
            } else if (key === 'Tab') {
                handled = true;
                const delta = isShift ? -1 : 1;
                newCol = selectedCol + delta;
                if (newCol < 0) {
                    if (newRow > 0) {
                        newRow -= 1;
                        newCol = totalCols - 1;
                    } else {
                        newCol = 0;
                    }
                } else if (newCol >= totalCols) {
                    if (newRow < totalRows - 1) {
                        newRow += 1;
                        newCol = 0;
                    } else {
                        newCol = totalCols - 1;
                    }
                }
            } else if (key === 'Enter') {
                handled = true;
                // Excel-style: move down
                newRow = Math.min(totalRows - 1, selectedRow + 1);
            }

            if (handled) {
                e.preventDefault();
                setSelection(newRow, newCol);
                ensureCellVisible();
            }
        });

        // Column resizing

        function attachResizeHandlers() {
            const handles = document.querySelectorAll('.col-resize-handle');
            const colgroup = document.getElementById('colgroup');
            let resizing = false;
            let startX = 0;
            let startWidth = 0;
            let targetColIndex = -1;

            function onMouseMove(e) {
                if (!resizing) return;
                const dx = e.clientX - startX;
                const newWidth = Math.max(40, startWidth + dx);
                const col = colgroup.children[targetColIndex];
                if (col) {
                    col.style.width = newWidth + 'px';
                }
            }

            function onMouseUp() {
                if (!resizing) return;
                resizing = false;
                document.removeEventListener('mousemove', onMouseMove);
                document.removeEventListener('mouseup', onMouseUp);
            }

            handles.forEach(handle => {
                handle.addEventListener('mousedown', (e) => {
                    e.preventDefault();
                    const idx = parseInt(handle.dataset.colIndex, 10);
                    const th = handle.parentElement;
                    const col = colgroup.children[idx];
                    if (!th || !col) return;

                    resizing = true;
                    targetColIndex = idx;
                    startX = e.clientX;
                    startWidth = th.getBoundingClientRect().width;

                    document.addEventListener('mousemove', onMouseMove);
                    document.addEventListener('mouseup', onMouseUp);
                });
            });
        }

        // Saving

        async function postJSON(url, obj) {
            const res = await fetch(url, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(obj)
            });
            return await res.json();
        }

        document.getElementById('saveBtn').addEventListener('click', async () => {
            // If currently editing, commit that edit before saving
            if (editing) {
                finishEdit(true);
            }
            const payload = { edits: edits };
            const res = await postJSON('/save', payload);
            if (res.ok) {
                alert('Changes saved.');
                edits = [];
            }
        });

        // Notify backend when tab closes (discarding unsaved edits)
        window.addEventListener("beforeunload", () => {
            navigator.sendBeacon("/close", JSON.stringify({}));
        });

        // Virtual scroll listener
        document.getElementById('grid').addEventListener('scroll', updateVisibleRows);

        loadData();
        </script>

        </body>
        </html>"""

    # ============================================================
    # CORE STRUCTURE AND METADATA
    # ============================================================

    def get_dtypes(self):
        """Return the data types of each column as a list."""
        return copy.deepcopy(self._dtypes)
    
    def set_dtypes(self,dtypes):
        """Set the data type for each column. Column values will be converted to the specified data type.
        
        Parameters:
        
        dtypes (list): a list of data type strings for each column. Can be specified as "integer" for integer data, "float" for float data, "string" for string data, or "bool" for Boolean data."""
        if len(dtypes) != self.num_columns():
            raise ValueError("Length of dtypes must match number of columns.")
        valid = ["integer", "float", "string", "bool"]
        for c, dtype in enumerate(dtypes):
            if dtype not in valid:
                raise ValueError(f"Invalid dtype '{dtype}'. Must be one of {valid}.")
            if dtype != self._dtypes[c]:
                self.convert_dtype(c, dtype)
        return self

    def detect_dtypes(self):
        """Automatically detect the data type of each column and set the Table data types."""
        self._dtypes.clear()
        for c in range(self.num_columns()):
            self._dtypes.append(self._return_col_type(self.get_column(c)))
        return self

    def convert_dtype(self,column,dtype):
        """Convert the data type of a single column to another data type.
        
        Parameters:

        column (integer, string): the index or header of the column to be converted.
        
        dtype (string): the data type to convert the column to. Can be specified as "integer" for integer data, "float" for float data, "string" for string data, or "bool" for Boolean data."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")

        if isinstance(column, str):
            column = self._headers.index(column)

        # Add bool support
        valid = {
            "integer": int,
            "float": float,
            "string": str,
            "bool": None  # custom handling below
        }

        if dtype not in valid:
            raise ValueError(f"Invalid dtype '{dtype}'. Must be one of {set(valid)}.")

        self._dtypes[column] = dtype

        # Custom boolean caster
        if dtype == "bool":
            true_set = self._boolean_true_values
            false_set = self._boolean_false_values

            def caster(val):
                # Already a bool
                if isinstance(val, bool):
                    return val

                # Numeric booleans
                if isinstance(val, int) and not isinstance(val, bool):
                    if val == 1:
                        return True
                    if val == 0:
                        return False

                # String booleans
                if isinstance(val, str):
                    low = val.strip().lower()
                    if low in true_set:
                        return True
                    if low in false_set:
                        return False

                # Unconvertible, leave unchanged
                return val

        else:
            # Use built-in caster for int, float, string
            caster = valid[dtype]

        # Apply conversion
        for row in self._data:
            try:
                row[column] = caster(row[column])
            except Exception:
                pass
        return self

    def autogenerate_headers(self):
        """Autogenerate generic headers for the data if there are no column headers."""
        if self.get_headers():
            return

        # Validate column count
        nc = self.num_columns()
        if nc == 0:
            raise ValueError("Cannot autogenerate headers for a table with no columns.")

        # Generate headers
        new_headers = [f"Column_{i+1}" for i in range(nc)]
        self.set_headers(new_headers=new_headers)
        return self

    def get_headers(self):
        """Return the columns headers."""
        return copy.deepcopy(self._headers)
    
    def set_headers(self,new_headers):
        """Insert new column headers ('new_headers') for all columns. If column headers already exist, this method replaces the existing headers. Headers may be defined for an empty Table.

        Parameters:

        new_headers (list): a list of the new headers for each column. If headers are not strings, they will be converted to strings."""
        # Validate type
        if not isinstance(new_headers, list):
            raise TypeError("Headers must be provided as a list.")
        
        if self.num_rows() == 0: #if this is an empty table
            self._headers = new_headers[:]
            return

        # Validate length
        if len(new_headers) != self.num_columns():
            raise ValueError(
                f"Header count ({len(new_headers)}) does not match number of columns ({self.num_columns()}).")

        # Convert all headers to strings
        self._headers = [str(h) for h in new_headers]
        return self
        
    def rename_column(self,column,new_header):
        """Replace an existing column header ('column') with a new header ('new_header').

        Parameters:

        column (integer, string): the index or header of the column to give a new header to.

        new_header (string): the new header for the column. If new_header is not a string it will be converted to a string."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in the table.")

        # Resolve column index
        if isinstance(column, str):
            column = self._headers.index(column)

        # Convert new header to string
        new_header = str(new_header)

        # Replace using a safe copy
        headers = self.get_headers()
        headers[column] = new_header

        # Delegate to set_headers for consistency
        self.set_headers(headers)
        return self
    
    def remove_header_space(self):
        """Replace spaces in header names with underscores. This is necessary when writing data to an SQL database table if the headers have spaces."""
        headers = self.get_headers()

        # Convert all headers to strings and replace spaces
        cleaned = [str(h).replace(" ", "_") for h in headers]

        self.set_headers(cleaned)
        return self

    def is_empty(self):
        """Returns True if this Table object has no data and False otherwise."""
        return not self._data

    def clear_data(self):
        """Remove the data, headers, and data types stored in the Table object."""
        self._data = []
        self._headers = []
        self._dtypes = []
        return self
    
    def copy(self):
        """Copy the headers, data types, and data of the Table and return a new Table object."""
        _out = TableTools().new_table()
        _out.set_data(self.get_data())
        _out.set_dtypes(self.get_dtypes())
        _out.set_headers(self.get_headers())
        return _out

    def num_rows(self):
        """Return the total number of rows in the data"""
        return len(self._data)

    def num_columns(self):
        """Return the total number of columns in the data."""
        if self._data:
            return len(self._data[0])
        else:
            return 0 

    def dimensions(self):
        """Return a tuple of (rows, columns), representing the total number of rows and columns in the data."""
        return (self.num_rows(),self.num_columns())

    # ============================================================
    # DATA ACCESS AND EXTRACTION
    # ============================================================

    def get_data(self):
        """Return all rows of data in the Table object as a nested list, not including column headers."""
        return copy.deepcopy(self._data)

    def set_data(self,rows):
        """Set the Table object data to a nested list of rows.
        
        Parameters: 
        
        rows (list): the list of rows."""
        self._data = copy.deepcopy(rows)
        self.detect_dtypes()
        return self

    def get_value(self,row,column):
        """Return a value at the specified position in the data.

        Parameters:

        row (integer): the index of the row of the value to be returned.

        column (integer. string): the index or header of the column of the value to be returned."""
        if row < 0 or row >= self.num_rows():
            raise ValueError("Specified row index is out of range.")

        # Validate column
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")

        # Resolve column index if header string
        if isinstance(column, str):
            column = self._headers.index(column)

        return self._data[row][column]
    
    def set_value(self,row,column,new_val):
        """Replace an existing value at a specific position in the data with a new value ('new_val').
        
        Parameters:

        row (integer): the index of the row of the value to be modified.

        column (integer, string): the index or header of the column of the value to be modified.

        new_val (integer,float,string): the new value to replace the existing value at the given position."""
        if row < 0 or row >= self.num_rows():
            raise ValueError("Specified row index is out of range.")

        # Validate column
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")

        # Resolve column index if header string
        if isinstance(column, str):
            column = self._headers.index(column)

        self._data[row][column] = new_val
        return self

    def get_row(self,row):
        """Return the values of a specified row as a list.

        Parameters:

        row (integer): the index of the row to be returned."""
        if row < 0 or row >= self.num_rows():
            raise ValueError("Specified index is out of range.")
        return self._data[row][:]

    def get_rows(self,rows):
        """Return a list of rows, each a list of values.

        Parameters:

        rows (list): a list of indices of the rows to be returned."""
        n_rows = self.num_rows()
        for r in rows:
            if r < 0 or r >= n_rows:
                raise ValueError("Specified index is out of range.")
            
        out_rows = []
        for r in rows:
            out_rows.append(self.get_row(r))
        return out_rows

    def get_column(self,column):
        """Return the values of a specified column as a list.

        Parameters:

        column (integer, string): the index or header of the column to be returned."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")

        if isinstance(column, str):
            column = self._headers.index(column)

        return [row[column] for row in self._data]
    
    def get_columns(self,columns):
        """Return the values of each specified column as a list.

        Parameters:

        column (list): a list of the indices or headers of the columns to be returned."""
        out_cols = []
        for c in columns:
            out_cols.append(self.get_column(c))
        return out_cols

    def value_in_row(self,row,value):
        """Determine if the specified value ('value') exists in the specified row ('row'). Returns True or False.
        
        Parameters:
        
        row (integer): the row to check for the specified value.
        
        value (integer, float, string): the value to check for."""
        if row < 0 or row >= self.num_rows():
            raise ValueError("Specified row index is out of range.")

        return value in self.get_row(row)

    def value_in_column(self,column,value):
        """Determine if the specified value ('value') exists in the specified column ('column'). Returns True or False.
        
        Parameters:
        
        column (integer, string): the index or header of the column to check for the specified value.
        
        value (integer, float, string): the value to check for."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")

        return value in self.get_column(column)

    def row_index_of_val(self,column,value):
        """Return the row index of a specified value ('value') in a specified column ('column'). Only the first index where the value appears will be returned.
        
        Parameters:
        
        column (integer, string): the index or header of the column the operation will be performed on.
        
        value (integer, float, string): the value to search for. Only the first index of where the value appears will be returned."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")
        if isinstance(column,str):
            column = self._headers.index(column)

        if not self.value_in_column(column,value):
            raise ValueError(f"Value '{value}' does not exist in this column.")
        for ix, row in enumerate(self._data):
            if row[column] == value:
                return ix
            
    def row_indices_of_val(self,column,value):
        """Return the row indices of a specified value ('value') in a specified column ('column'). All row indices where the value appears will be returned in a list.
        
        Parameters:
        
        column (integer, string): the index or header of the column the operation will be performed on.
        
        value (integer, float, string): the value to search for. All indices of the target value will be returned."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")
        if isinstance(column,str):
            column = self._headers.index(column)

        if not self.value_in_column(column,value):
            raise ValueError(f"Value '{value}' does not exist in this column.")
        
        indices = []
        for ix, row in enumerate(self._data):
            if row[column] == value:
                indices.append(ix)
        return indices

    def column_index_of_val(self,row,value):
        """Return the column index of a specified value ('value') in a specified row ('row'). Only the first index where the value appears will be returned.
        
        Parameters:
        
        row (integer): the row to check for the specified value index.
        
        value (integer, float, string): the value to search for. Only the first index of where the value appears will be returned."""
        if row < 0 or row >= self.num_rows():
            raise ValueError("Specified row index is out of range.")

        if not self.value_in_row(row, value):
            raise ValueError(f"Value '{value}' does not exist in this row.")

        r = self.get_row(row)
        for c, val in enumerate(r):
            if val == value:
                return c
            
    def column_indices_of_val(self,row,value):
        """Return the column indices of a specified value ('value') in a specified row ('row'). All indices where the value appears will be returned as a list.
        
        Parameters:
        
        row (integer): the row to check for the specified value index.
        
        value (integer, float, string): the value to search for. All indices of the target value will be returned."""
        if row < 0 or row >= self.num_rows():
            raise ValueError("Specified row index is out of range.")

        if not self.value_in_row(row, value):
            raise ValueError(f"Value '{value}' does not exist in this row.")

        indices = []
        r = self.get_row(row)
        for c, val in enumerate(r):
            if val == value:
                indices.append(c)
        return indices

    def unique_vals_in_row(self,row):
        """Return a list of unique values in the specified row.
        
        Parameters:
        
        row (integer): the row the operation will be performed on."""
        if row < 0 or row >= self.num_rows():
            raise ValueError("Specified row index is out of range.")

        return list(set(self.get_row(row)))

    def get_sample(self,num_rows):
        """Get a random sample of rows ('num_rows') and return a new Table object. Column headers and data types are transferred to the new Table.
        
        Parameters:
        
        num_rows (integer): the number of rows in the new Table."""
        if num_rows > self.num_rows():
            num_rows = self.num_rows()
        if num_rows < 0:
            raise ValueError("num_rows must be non-negative.")

        _out = TableTools().new_table()

        samp = _Stats().get_random_sample(self.get_data(),num_rows)
        _out.set_data(samp)
        _out.set_headers(self.get_headers())
        _out.set_dtypes(self.get_dtypes())
        return _out

    def count_in_row(self,row,value):
        """Count and return the number of occurrences of a specified value ('value') in the specified row ('row').
        
        Parameters:
        
        row (integer): the row to count the specified value in.
        
        value (integer, float, string): the value to count."""
        if row < 0 or row >= self.num_rows():
            raise ValueError("Specified row index is out of range.")
        
        return self.get_row(row).count(value)
        
    def count_in_column(self,column,value):
        """Count and return the number of occurrences of a specified value ('value') in the specified column ('column').
        
        Parameters:
        
        column (integer, string): the index or header of the column to count the specified value in.
        
        value (integer, float, string): the value to count."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")

        return self.get_column(column).count(value)

    def unique_vals_in_col(self,column):
        """Return a list of unique values in the specified column.
        
        Parameters:
        
        column (integer, string): the index or header of the column the operation will be performed on."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")

        return list(set(self.get_column(column)))

    # ============================================================
    # DISPLAY AND INSPECTION
    # ============================================================

    def display_object_functions(self):
        """Print a list of all Table object functions."""
        print(f"Number of {__class__.__name__[1:]} functions: {len([f for f in dir(__class__) if not f.startswith('__') and not f.startswith('_')])}")
        for f in [f for f in dir(__class__) if not f.startswith("__") and not f.startswith('_')]:
            print(f)

    def display_rows(self,start_row=None,end_row=None):
        """Print a range of rows. The specified range is inclusive.

        Parameters:

        start_row (integer): the index of the first row to be displayed. If None start_row defaults to 0.

        end_row (integer): the index of the last row to be displayed. The range specified by start_row and end_row is inclusive. If None end_row defaults to the length of the data."""
        if start_row is None:
            start_row = 0
        if end_row is None:
            end_row = self.num_rows() - 1

        
        if start_row < 0 or end_row < 0:
            raise ValueError("Row indices must be non-negative.")

        if start_row > end_row:
            raise ValueError("start_row cannot be greater than end_row.")

        if end_row >= self.num_rows():
            raise ValueError("Row index out of range.")

        nc = self.num_columns()
        if not self.get_headers():
            self.autogenerate_headers()
        tmp = TableTools().new_table()
        tmp.set_data(self.get_data()[start_row:end_row+1])
        tmp.set_headers(self.get_headers())
        to_print = []
        for c in range(nc):
            col = tmp.get_column(c) #get each column as a list
            col.insert(0,tmp._headers[c]) #insert the header at the start of the list
            col = [str(v) for v in col] #convert each value to a string
            lvs = [len(v) for v in col] #find the length of each value, including the header
            ml = max(lvs) #find the max of the lengths of this column
            row = [col[i]+(" " * (ml-lvs[i]+2)) for i in range(len(col))] #add a number of blank spaces to the end of each value to equal the max column length + 2, so the columns are aligned
            row.insert(1,'-'*(ml+2)) #insert a row of dashes to separate the header from the values
            to_print.append(row)
        for r in range(len(to_print[0])):
            row = []
            for c in range(nc):
                row.append(to_print[c][r]) #convert the corresponding values of each column into a row
            row = ''.join(row)
            print(row)
        print("")

    def display_columns(self,start_col=None,end_col=None):
        """Print a range of columns. The specified range is inclusive.

        Parameters:

        start_col (integer): the index of the first column to be displayed. If None start_col defaults to 0.

        end_col (integer): the index of the last column to be displayed. The range specified by start_col and end_col is inclusive. If None end_col defaults to the number of the columns."""
        # Default to full column range
        if start_col is None:
            start_col = 0
        if end_col is None:
            end_col = self.num_columns() - 1

        # Validate indices
        if not isinstance(start_col, int) or not isinstance(end_col, int):
            raise TypeError("Column indices must be integers.")

        if start_col < 0 or end_col < 0:
            raise ValueError("Column indices must be non-negative.")

        if start_col > end_col:
            raise ValueError("start_col cannot be greater than end_col.")

        if end_col >= self.num_columns():
            raise ValueError("Column index out of range.")

        # Ensure headers exist
        if not self.get_headers():
            self.autogenerate_headers()

        nr = self.num_rows()

        # Build a temporary table
        tmp = TableTools().new_table()
        tmp.set_data(self.get_data())
        tmp.set_headers(self.get_headers())

        to_print = []

        # Build each column block
        for c in range(start_col, end_col + 1):
            col = tmp.get_column(c)
            col = [tmp._headers[c]] + col
            col = [str(v) for v in col]
            lengths = [len(v) for v in col]
            max_len = max(lengths)
            padded = [
                col[i] + " " * (max_len - lengths[i] + 2)
                for i in range(len(col))
            ]
            padded.insert(1, "-" * (max_len + 2))
            to_print.append(padded)

        if not to_print:
            print("")
            return

        # Print row by row (header + separator + all data rows)
        for r in range(len(to_print[0])):
            row = ''.join(to_print[c][r] for c in range(len(to_print)))
            print(row)

        print("")

    def display_head(self, num_rows=10):
        """Print the first rows ('num_rows') of the Table.

        Parameters:

        num_rows (integer): the number of first rows to be printed."""
        end_row = min(num_rows - 1, self.num_rows() - 1)
        self.display_rows(0, end_row)
    
    def display_tail(self, num_rows=10):
        """Print the last rows ('num_rows') of the Table.

        Parameters:

        num_rows (integer): the number of last rows to be printed."""
        last = self.num_rows() - 1
        start = max(0, last - num_rows + 1)
        self.display_rows(start, last)

    def display_sample(self, num_rows=10):
        """Extract and print a sample of rows ('num_rows').

        Parameters:

        num_rows (integer): the number of first rows to be printed."""
        sample = self.get_sample(num_rows)
        sample.display_rows()

    def info(self):
        """Print a summary of the Table object, including number of rows, number of columns, column headers, data types, and memory footprint."""
        import sys

        def _get_size(obj, seen=None):
            """Recursively calculate memory footprint in bytes."""
            size = sys.getsizeof(obj)
            if seen is None:
                seen = set()
            obj_id = id(obj)
            if obj_id in seen:
                return 0
            seen.add(obj_id)

            if isinstance(obj, dict):
                size += sum((_get_size(v, seen) for v in obj.values()))
                size += sum((_get_size(k, seen) for k in obj.keys()))
            elif hasattr(obj, '__dict__'):
                size += _get_size(obj.__dict__, seen)
            elif isinstance(obj, (list, tuple, set, frozenset)):
                size += sum((_get_size(i, seen) for i in obj))
            return size

        # calculate memory footprint of the table’s core data
        mem_bytes = (
        _get_size(self._data) +
        _get_size(self._headers) +
        _get_size(self._dtypes))
        mem_kb = mem_bytes / 1024
        mem_mb = mem_kb / 1024

        print("\n" + "="*40)
        print("Table Info")
        print("="*40)
        print(f"Rows            : {self.num_rows()}")
        print(f"Columns         : {self.num_columns()}")
        print(f"Headers         : {', '.join(self.get_headers())}")
        print(f"Data Types      : {', '.join(self.get_dtypes())}")
        print(f"Memory Footprint: {mem_bytes:,} bytes "
            f"({mem_kb:.2f} KB / {mem_mb:.2f} MB)")
        print("="*40 + "\n")

    def open_viewer(self):
        """Open an interactive browser-based data viewer with spreadsheet-style navigation and cell editing. Script execution pauses until the browser is closed. 
        Edits are not saved to the Table data unless explicitly saved by clicking the Save button. Closing the browser discards any unsaved edits.
        
        Navigation:
        - Click once to select a cell
        - Double-click (or press F2) to edit the selected cell
        - Arrow keys move between cells (when not editing)
        - Enter moves down; Shift+Enter moves up
        - Tab moves right; Shift+Tab moves left
        - While editing, arrow keys move the text cursor within the cell
        - Press Enter to commit an edit, or Esc to cancel it
        - Drag the small handle on the right edge of any column header to resize it."""
        import json
        import threading
        import socket
        import webbrowser
        from http.server import HTTPServer, BaseHTTPRequestHandler

        state = {
            "closed": threading.Event(),
        }

        def find_free_port(start=8765, end=8900):
            for port in range(start, end + 1):
                with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                    try:
                        s.bind(("127.0.0.1", port))
                        return port
                    except OSError:
                        continue
            raise RuntimeError("No free port available for viewer.")

        port = find_free_port()
        table_ref = self

        class ViewerHandler(BaseHTTPRequestHandler):

            def log_message(self, format, *args):
                return

            def _send_json(self, obj):
                payload = json.dumps(obj).encode("utf-8")
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(payload)))
                self.end_headers()
                self.wfile.write(payload)

            def _send_html(self, html):
                body = html.encode("utf-8")
                self.send_response(200)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)

            def do_GET(self):
                if self.path == "/":
                    self._send_html(table_ref._viewer_html())

                elif self.path == "/data":
                    self._send_json({
                        "headers": table_ref._headers,
                        "rows": table_ref._data
                    })

                else:
                    self.send_error(404)

            def do_POST(self):
                length = int(self.headers.get("Content-Length", "0"))
                raw = self.rfile.read(length)
                body = json.loads(raw.decode("utf-8")) if raw else {}

                if self.path == "/save":
                    edits = body.get("edits", [])

                    # Apply edits immediately and cumulatively
                    for edit in edits:
                        r = edit["row"]
                        c = edit["col"]
                        v = edit["value"]
                        table_ref._data[r][c] = v

                    self._send_json({"ok": True})

                elif self.path == "/close":
                    state["closed"].set()
                    self._send_json({"ok": True})

                    # Clean shutdown from another thread
                    threading.Thread(target=server.shutdown, daemon=True).start()

                else:
                    self.send_error(404)

        server = HTTPServer(("127.0.0.1", port), ViewerHandler)

        def serve():
            server.serve_forever()
            server.server_close()

        threading.Thread(target=serve, daemon=True).start()

        webbrowser.open(f"http://127.0.0.1:{port}/")

        # Wait until the browser/tab signals closure
        state["closed"].wait()

        self.detect_dtypes()

    # ============================================================
    # COLUMN OPERATIONS
    # ============================================================

    def append_column(self,vals,col_header):
        """Append a new column of values ('vals') to the end of the data. This function can also be used to add a first column of data to a new Table object.

        Parameters:

        vals (list): the column of values to be appended.
        
        col_header (string): the header of the column to be appended."""
        ld = len(self._data)
        lv = len(vals)

        if ld > 0 and lv != ld:
            raise ValueError("New column must be the same length as existing rows.")

        v_vals = vals[:]
        dtype = self._return_col_type(vals)

        if not self._headers and not self.is_empty():
            self.autogenerate_headers()

        self._headers.append(col_header)
        self._dtypes.append(dtype)

        if self._data:  # table already has rows
            for i in range(ld):
                self._data[i].append(v_vals[i])
        else:  # empty table, build new rows
            for val in v_vals:
                self._data.append([val])
        return self

    def append_columns(self,vals,col_headers):
        """Append a new column of values ('vals') to the end of the data. This function can also be used to add a first column of data to a new Table object.

        Parameters:

        vals (list): a nested list of column values to be appended.
        
        col_headers (list): a list of headers of the columns to be appended."""
        if len(vals) != len(col_headers):
            raise ValueError("Number of value lists must match number of headers.")

        for v, h in zip(vals, col_headers):
            self.append_column(v, h)
        return self
 
    def insert_column(self,vals,index,col_header):
        """Insert a new column of values ('vals') in the data at the specified index ('index'). This function can also be used to add a first column of data to a new Table object.

        Parameters:

        vals (list): the column of values to be inserted.

        index (integer): the index where the column will be inserted. Columns after the inserted column will be shifted to the right.

        col_header (string): the name of the header for the new column."""
        ld = len(self._data)
        lv = len(vals)

        if ld > 0 and lv != ld:
            raise ValueError("New column must be the same length as existing rows.")

        v_vals = vals[:]
        dtype = self._return_col_type(vals)

        if not self._headers and not self.is_empty():
            self.autogenerate_headers()

        self._headers.insert(index, col_header)
        self._dtypes.insert(index, dtype)

        if self._data:  # table already has rows
            for i in range(ld):
                self._data[i].insert(index, v_vals[i])
        else:  # empty table, build new rows
            for val in v_vals:
                self._data.append([val])
        return self

    def insert_columns(self, vals, index, col_headers):
        """Insert new columns of values ('vals') in the data at the specified index ('index'). This function can also be used to add the first columns of data to a new Table object.

        Parameters:

        vals (list): a nested list of column values to be inserted. Each inner list must match the number of rows if the table is not empty.

        index (integer): the index where the first column will be inserted. Columns after the inserted columns will be shifted to the right.

        col_headers (list): a list of headers for the new columns."""
        ld = len(self._data)

        if len(vals) != len(col_headers):
            raise ValueError("Number of value lists must match number of headers.")

        for v in vals:
            if ld > 0 and len(v) != ld:
                raise ValueError("Each new column must be the same length as existing rows.")

        if not self._headers and not self.is_empty():
            self.autogenerate_headers()

        # Insert headers and dtypes in order
        for i, (v, h) in enumerate(zip(vals, col_headers)):
            dtype = self._return_col_type(v)
            self._headers.insert(index + i, h)
            self._dtypes.insert(index + i, dtype)

        if self._data:  # table already has rows
            for row_idx in range(ld):
                for i, v in enumerate(vals):
                    self._data[row_idx].insert(index + i, v[row_idx])
        else:  # empty table, build new rows
            for row_vals in zip(*vals):
                self._data.append(list(row_vals))
        return self

    def duplicate_column(self,column,duplicates):
        """Duplicate a specified column of data.

        Parameters:

        column (integer, string): the index or header of the column to be duplicated. Columns after the duplicated columns will be shifted to the right.

        duplicates (integer): the number of times to duplicate the specified column."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")

        if isinstance(column, str):
            column = self._headers.index(column)

        if not self._headers:
            self.autogenerate_headers()

        header = self._headers[column]
        dtype = self._dtypes[column]

        # Duplicate headers and dtypes
        for _ in range(duplicates):
            self._headers.insert(column, header)
            self._dtypes.insert(column, dtype)

        # Duplicate column values in each row
        for row in self._data:
            val = row[column]
            for _ in range(duplicates):
                row.insert(column, val)
        return self

    def remove_columns(self,start_col,end_col):
        """Remove a range of columns from the data.

        Parameters:

        start_col (integer): the index of the first column to be removed.

        end_col (integer): the index of the last column to be removed. The range specified by start_col and end_col is inclusive."""
        if not isinstance(start_col, int) or not isinstance(end_col, int):
            raise TypeError("Column indices must be integers.")

        if start_col < 0 or end_col < 0:
            raise ValueError("Column indices must be non-negative.")

        if start_col > end_col:
            raise ValueError("start_col cannot be greater than end_col.")

        if end_col >= self.num_columns():
            raise ValueError("Column index out of range.")

        # Update headers and dtypes
        self._headers = self._headers[:start_col] + self._headers[end_col+1:]
        self._dtypes = self._dtypes[:start_col] + self._dtypes[end_col+1:]

        # Remove the specified columns from each row
        for i in range(len(self._data)):
            self._data[i] = self._data[i][:start_col] + self._data[i][end_col+1:]
        return self

    def remove_duplicate_columns(self):
        """Remove all duplicate columns from the data."""
        headers = self.get_headers()
        dtypes = self.get_dtypes()
        nc = self.num_columns()
        cols = [self.get_column(c) for c in range(nc)]

        unique_cols = []
        unique_headers = []
        unique_dtypes = []

        for c in range(nc):
            if cols[c] not in unique_cols:
                unique_cols.append(cols[c])
                unique_headers.append(headers[c])
                unique_dtypes.append(dtypes[c])

        # Rebuild data with only unique columns
        self._data = [list(row) for row in zip(*unique_cols)]
        self._headers = unique_headers
        self._dtypes = unique_dtypes
        return self

    def drop_column(self, column):
        """Remove a single specified column from the data.

        Parameters:

        column (integer, string): the index or header of the column to be removed."""

        # Validate column existence
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")

        # Convert header to index
        if isinstance(column, str):
            column = self._headers.index(column)

        # Delegate to existing range-based remover
        self.remove_columns(column, column)
        return self

    def shift_column(self,column,new_index):
        """Move a specified column from its current position in the data to a new position as specified by an index ('new_index').
        
        Parameters:
        
        column (integer, string): the index or header of the column to be moved.
        
        new_index (integer): the index of the new position for the specified column."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")

        if isinstance(column, str):
            column = self._headers.index(column)

        header = self._headers[column]
        c = self.get_column(column)

        self.remove_columns(column, column)
        self.insert_column(c, new_index, header)
        return self

    def shift_columns(self, columns, new_index):
        """Move multiple specified columns from their current positions in the data to a new position as specified by an index ('new_index').
        
        Parameters:
        
        columns (list): a list of indices of the columns to be moved.
        
        new_index (integer): the index of the new position for the columns. If equal to the number of columns, the columns will be moved to the end."""
        n_cols = self.num_columns()
        for c in columns:
            if c < 0 or c >= n_cols:
                raise ValueError("Specified column index is out of range.")
        if new_index < 0 or new_index > n_cols:
            raise ValueError("Specified index 'new_index' is out of range.")

        # Extract the columns to move
        moved = [self.get_column(c) for c in columns]
        moved_headers = [self._headers[c] for c in columns]

        # Remove them from the table (sorted descending so indices stay valid)
        for c in sorted(columns, reverse=True):
            self.remove_columns(c, c)

        # Insert them back starting at new_index
        for i, (col_vals, header) in enumerate(zip(moved, moved_headers)):
            self.insert_column(col_vals, new_index + i, header)
        return self

    def replace_column(self,column,vals,col_header):
        """Replace an existing column ('column') with a new column ('vals').

        Parameters:

        column (integer, string): the index or header of the column to be replaced.

        vals (list): the list of values to replace the old column.

        col_header (string): the name of the header for the new column."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")

        ld = len(self._data)
        lv = len(vals)
        if lv != ld:
            raise ValueError("New column must be the same length as existing rows.")

        if isinstance(column, str):
            column = self._headers.index(column)

        v_vals = vals[:]
        for i in range(ld):
            self._data[i][column] = v_vals[i]

        self._headers[column] = col_header
        self._dtypes[column] = self._return_col_type(vals)
        return self
    
    def pop_column(self,column):
        """Return the values of the specified column as a list and remove the column from the data.
        
        Parameters:

        column (integer, string): the index or header of the column to be returned."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")

        if isinstance(column, str):
            column = self._headers.index(column)

        col = self.get_column(column)
        self.remove_columns(column, column)
        return col

    def swap_columns(self, col1, col2):
        """Swap the positions of two specified columns in the data.

        Parameters:

        col1 (integer, string): the index or header of the first column to be swapped.

        col2 (integer, string): the index or header of the second column to be swapped."""
        if not self.column_exists(col1) or not self.column_exists(col2):
            raise ValueError("One or both specified columns do not exist in Table.")

        if isinstance(col1, str):
            col1 = self._headers.index(col1)
        if isinstance(col2, str):
            col2 = self._headers.index(col2)

        # Swap headers and dtypes
        self._headers[col1], self._headers[col2] = self._headers[col2], self._headers[col1]
        self._dtypes[col1], self._dtypes[col2] = self._dtypes[col2], self._dtypes[col1]

        # Swap values in each row
        for row in self._data:
            row[col1], row[col2] = row[col2], row[col1]
        return self

    def reverse_columns(self):
        """Reverse the order of all columns in the data."""
        self._headers.reverse()
        self._dtypes.reverse()
        for i in range(len(self._data)):
            self._data[i].reverse()
        return self

    def sort_columns(self, method):
        """Sort all columns in the Table by their headers using the specified method ('method') and return a new Table.
        
        Parameters:
        
        method (string): the method to sort by. Can be specified as "ascending", "asc", "low" or "descending", "desc", "high"."""
        method = method.lower()

        ascending_aliases = ("ascending", "asc", "low")
        descending_aliases = ("descending", "desc", "high")

        if method in ascending_aliases:
            reverse = False
        elif method in descending_aliases:
            reverse = True
        else:
            raise ValueError(
                "Method must be one of: 'ascending', 'descending', 'asc', 'desc', 'low', or 'high'."
            )

        # Pair headers, dtypes, and column values together
        cols = list(zip(self._headers, self._dtypes, zip(*self._data)))

        # Sort by header
        cols.sort(key=lambda x: x[0], reverse=reverse)

        # Build output table
        _out = TableTools().new_table()
        _out.set_data([list(row) for row in zip(*[c for _, _, c in cols])])
        _out.set_headers([h for h, _, _ in cols])
        _out.set_dtypes([d for _, d, _ in cols])
        return _out

    def column_exists(self,column):
        """Determine if the specified column ('column') exists in the Table object. Returns True or False.

        Parameters:
        
        column (integer, string): the index or header of the column to be checked."""
        if isinstance(column, str):
            return column in self._headers

        # Integer index lookup
        if isinstance(column, int):
            if column < 0:
                return False
            return column < self.num_columns()

        # Anything else, invalid
        return False
    
    def column_index(self,column):
        """Return the index of a specified column ('column') if it exists in the Table object.
        
        Parameters:
        
        column (string): the header of the column to determine the index of."""
        if not isinstance(column, str):
            raise TypeError("column_index() requires a column header (string).")

        if column not in self._headers:
            raise ValueError(f"Column '{column}' does not exist in Table.")

        return self._headers.index(column)

    def column_subset(self,columns):
        """Return a new Table object that includes only the specified columns ('columns').
        
        Parameters:
        
        columns (list): the list of column headers or columns indices to subset in the output Table."""
        if not isinstance(columns, list):
            raise TypeError("columns must be a list of column headers or indices.")

        # Normalize all columns to indices
        indices = []
        for col in columns:
            if not self.column_exists(col):
                raise ValueError(f"Column '{col}' does not exist in Table.")

            if isinstance(col, int):
                indices.append(col)
            else:  # string
                indices.append(self.column_index(col))

        # Build new table
        _out = TableTools().new_table()

        for idx in indices:
            col_data = self.get_column(idx)
            header = self._headers[idx]
            _out.append_column(col_data, header)

        return _out

    # ============================================================
    # ROW OPERATIONS
    # ============================================================

    def append_row(self,vals):
        """Append a new row of values ('vals') to the end of the data. Can be used to add a first row to an empty Table.

        Parameters:

        vals (list): the row of values to be appended. If the Table is not empty, the new row must be the same length as the number columns."""
        if not isinstance(vals,list):
            raise ValueError("'vals' must be a list of values.")
        nc = self.num_columns()
        if len(vals) != nc and nc > 0: #the check for zero columns is useful if appending to an empty table.
            raise ValueError("New row must be the same length as existing rows.")

        v = vals[:]
        if self.is_empty():
            self._row_dtypes(v)
        self._data.append(v)
        return self
    
    def append_rows(self,vals):
        """Append a list of rows ('vals') to the end of the data.

        Parameters:

        vals (list): the list of rows to be appended. If the Table is not empty, each row must be the same length as the number columns."""
        if not isinstance(vals,list):
            raise ValueError("'vals' must be a list of values.")
        for v in vals:
            self.append_row(v)
        return self

    def insert_row(self,vals,index):
        """Insert a new row of values ('vals') in the data at the specified index ('index').

        Parameters:

        vals (list): the row of values to be inserted. If the Table is not empty, new row must be the same length as the number of columns.

        index (integer): the index where the row will be inserted. Rows after the inserted row will be shifted down."""
        if not isinstance(vals,list):
            raise ValueError("'vals' must be a list of values.")
        nc = self.num_columns()
        if len(vals) != nc and nc > 0: #the check for zero columns is useful if appending to an empty table.
            raise ValueError("New row must be the same length as existing rows.")
        if index < 0 or index > self.num_rows():
            raise ValueError("Specified index is out of range.")
        
        v = vals[:]
        if self.is_empty():
            self._row_dtypes(v)
        self._data.insert(index,v)
        return self
    
    def insert_rows(self, vals, index):
        """Insert multiple rows of values ('vals') in the data at the specified index ('index').
        
        Parameters:
        
        vals (list): the list of rows to be inserted. If the Table is not empty, each row must be the same length as the number of columns.
        
        index (integer): the index where the rows will be inserted. Rows at and after the inserted rows will be shifted down.
        """
        for v in vals:
            self.insert_row(v, index)
            index += 1
        return self

    def duplicate_row(self,row,duplicates):
        """Duplicate a specified row of data.

        Parameters:

        row (integer): the index of the row to be duplicated. Rows after the duplicated rows will be shifted down.

        duplicates (integer): the number of times to duplicate the specified row."""
        if row < 0 or row >= self.num_rows():
            raise ValueError("Specified index is out of range.")
        for d in range(duplicates):
            self._data.insert(row+d+1,self._data[row][:])
        return self

    def remove_rows(self,start_row,end_row):
        """Remove a range of rows from the data.

        Parameters:

        start_row (integer): the index of the first row to be removed.

        end_row (integer): the index of the last row to be removed. The range specified by start_row and end_row is inclusive."""
        if start_row < 0 or end_row < 0 or start_row >= self.num_rows() or end_row >= self.num_rows():
            raise ValueError("Specified row indices are out of range.")

        if start_row > end_row:
            raise ValueError("start_row cannot be greater than end_row.")
        
        self._data = self._data[:start_row] + self._data[end_row+1:]
        return self
 
    def remove_duplicate_rows(self):
        """Remove all duplicate rows from the data. The first occurrence of each duplicate is presereved."""
        no_dups = _ListOps().remove_duplicates_ordered(self._data[:])
        self.set_data(no_dups)
        return self

    def shift_row(self,row,new_index):
        """Move a specified row from its current position in the data to a new position as specified by an index ('new_index').
        
        Parameters:
        
        row (integer): the index of the row to be moved.
        
        new_index (integer): the index of the new position for the row."""
        if row == new_index:
            return
        if row < 0 or row >= self.num_rows():
            raise ValueError("Specified index is out of range.")
        if new_index < 0 or new_index >= self.num_rows():
            raise ValueError("Specified index 'new_index' is out of range.")
        r = self.get_row(row)
        self.remove_rows(row,row)
        self.insert_row(r,new_index)
        return self

    def shift_rows(self, rows, new_index):
        """Move multiple specified rows from their current positions in the data to a new position as specified by an index ('new_index').
        
        Parameters:
        
        rows (list): a list of indices of the rows to be moved.
        
        new_index (integer): the index of the new position for the rows. If equal to the number of rows,the rows will be moved to the end."""
        n_rows = self.num_rows()
        for r in rows:
            if r < 0 or r >= n_rows:
                raise ValueError("Specified row index is out of range.")
        if new_index < 0 or new_index > n_rows:
            raise ValueError("Specified index 'new_index' is out of range.")

        # Extract the rows to move
        moved = [self.get_row(r) for r in rows]

        # Remove them from the table (sorted descending so indices stay valid)
        for r in sorted(rows, reverse=True):
            self.remove_rows(r, r)

        # Insert them back starting at new_index
        for i, row in enumerate(moved):
            self.insert_row(row, new_index + i)
        return self

    def replace_row(self,row,vals):
        """Replace an existing row ('row') with a new row ('vals').

        Parameters:

        row (integer): the index of the row to be replaced.

        vals (list): the row of values to replace the old row. Row must match the number of columns."""
        if row < 0 or row >= self.num_rows():
            raise ValueError("Specified index is out of range.")
        if not isinstance(vals,list):
            raise ValueError("'vals' must be a list of values.")
        nc = self.num_columns()
        if len(vals) != nc and nc > 0: #the check for zero columns is useful if appending to an empty table.
            raise ValueError("New row must be the same length as existing rows.")
        
        v = vals[:]
        self._data[row] = v
        if self.num_rows() == 1:
            self.detect_dtypes()
        return self

    def pop_row(self,row):
        """Return the values of a specified row as a list and remove the row from the data.

        Parameters:

        row (integer): the index of the row to be returned."""
        if row < 0 or row >= self.num_rows():
            raise ValueError("Specified index is out of range.")
        r = self.get_row(row)
        self.remove_rows(row,row)
        return r

    def swap_rows(self, row1, row2):
        """Swap the positions of two specified rows in the data.
        
        Parameters:
        
        row1 (integer): the index of the first row to be swapped.
        
        row2 (integer): the index of the second row to be swapped."""
        n_rows = self.num_rows()
        if row1 < 0 or row1 >= n_rows:
            raise ValueError("Specified index 'row1' is out of range.")
        if row2 < 0 or row2 >= n_rows:
            raise ValueError("Specified index 'row2' is out of range.")
        if row1 == row2:
            return

        self._data[row1], self._data[row2] = self._data[row2], self._data[row1]
        return self

    def reverse_rows(self):
        """Reverse the order of all rows in the data."""
        self._data.reverse()
        return self

    def sort_rows(self,column,method):
        """Sort all rows in the Table by a specified column ('column') and method ('method') and return a new Table.
        
        Parameters:
        
        column (integer, string): the index or header of the column to sort columns by.
        
        method (string): the method to sort by. Can be specified as "ascending", "asc", "low" or "descending", "desc", "high"."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")

        # Convert header to index
        if isinstance(column, str):
            column = self._headers.index(column)

        # Normalize and validate method
        method = method.lower()

        ascending_aliases = ("ascending", "asc", "low")
        descending_aliases = ("descending", "desc", "high")

        if method in ascending_aliases:
            reverse = False
        elif method in descending_aliases:
            reverse = True
        else:
            raise ValueError(
                "Method must be one of: 'ascending', 'descending', 'asc', 'desc', 'low', or 'high'."
            )

        # Sort data
        data = self.get_data()  # safe copy
        try:
            data.sort(key=lambda row: row[column], reverse=reverse)
        except Exception as e:
            raise TypeError(f"Column '{column}' contains non-sortable values: {e}")

        # Build output table
        _out = TableTools().new_table()
        _out.set_data(data)
        _out.set_headers(self.get_headers())
        _out.set_dtypes(self.get_dtypes())
        return _out

    def dissolve_rows(self,dissolve_on,dissolve_columns = None,methods = None,sep = "", decimals = 2):
        """Aggregate rows with duplicate values in a specified column ('dissolve_on') into single rows and return a new Table. Values from other columns ('dissolve_columns') in the dissolved rows may be aggregated according to a specified method.
        
        Parameters:
        
        dissolve_on (integer, string): the index or header of the column with duplicate values to dissolve rows by.
        
        dissolve_columns (integer, string, list): the other column or columns in each duplicate row whose values will be combined according to a method. Can be left blank to include all columns. For a single column, can be specified as the index or header of a column. For multiple columns, can be specified as a list of column indices or headers. Columns not specified will be aggregated by the default "First" method.
        
        methods (string, list): If a single dissolve column is specified or no dissolve columns are specified (all columns will be included), this should be a string representing the method used to aggregate values. If multiple dissolve columns are specified, this should be a list of strings representing the aggregation method for each column. Possible aggregation methods are "First", "Last", "Count", "Unique", "Sum", "Min", "Max", "Range", "Mean", "Median", "Mode, "IQ1", "IQ3", "IQR", "Standard Deviation", "Standard Error", "Variance", "Coefficient of Variation", "Skewness", "Kurtosis", or "Concatenate". May be left blank to use "First" as the default aggregation method. The 'dissolve_on' column will always use the "First" method.
        
        sep (string): If the Concatenate aggregation method is used for any dissolve columns, a separator may be specified so that aggregated values can be retrieved as a string and split at the separator.
        
        decimals (integer): the number of decimal places float values will be rounded to."""
        if isinstance(dissolve_on, str):
            dissolve_on = self._headers.index(dissolve_on)

        num_cols = self.num_columns()
        headers = self.get_headers()
        dtypes = self.get_dtypes()
        data = self._data

        # Build groups: key -> list of rows
        groups = {}
        for row in data:
            key = row[dissolve_on]
            groups.setdefault(key, []).append(row)

        # Internal (non-stats) aggregation functions
        def _concatenate(col):
            return sep.join(str(v) for v in col)

        INTERNAL_AGG_FUNCS = {
            "First": lambda col: col[0],
            "Last": lambda col: col[-1],
            "Unique": lambda col: len(set(col)),
            "Concatenate": _concatenate,
        }
        
        # Stats-based aggregation methods (delegated to _Stats)
        STATS_METHODS = {
            "Count",
            "Unique Values",
            "First", 
            "Last",
            "Sum",
            "Max",
            "Min",
            "Range",
            "Mean",
            "Median",
            "Mode",
            "IQ1",
            "IQ3",
            "IQR",
            "Standard Deviation",
            "Standard Error",
            "Variance",
            "Coefficient of Variation",
            "Skewness",
            "Kurtosis",
        }

        # Normalize dissolve_columns and methods into per-column rules
        col_methods = {c: "First" for c in range(num_cols)}

        if dissolve_columns is None:
            if methods is not None:
                if not isinstance(methods, str):
                    raise ValueError("When dissolve_columns is None, 'methods' must be a single string or None.")
                for c in range(num_cols):
                    col_methods[c] = methods
        else:
            if isinstance(dissolve_columns, (int, str)):
                dissolve_columns = [dissolve_columns]

            # Convert headers to indices
            norm_cols = []
            for mc in dissolve_columns:
                if isinstance(mc, str):
                    norm_cols.append(self._headers.index(mc))
                else:
                    norm_cols.append(mc)
            dissolve_columns = norm_cols

            # Normalize methods
            if methods is None:
                for c in dissolve_columns:
                    col_methods[c] = "First"
            elif isinstance(methods, str):
                for c in dissolve_columns:
                    col_methods[c] = methods
            else:
                if len(methods) != len(dissolve_columns):
                    raise ValueError("Length of 'methods' must match length of 'dissolve_columns'.")
                for c, m in zip(dissolve_columns, methods):
                    col_methods[c] = m

        # Precompute aggregated rows for keys with duplicates
        agg_by_key = {}
        math_tool = _Math()
        stats_tool = _Stats()

        for key, rows in groups.items():
            if len(rows) == 1:
                continue

            agg_row = []
            for c in range(num_cols):
                col_vals = [r[c] for r in rows]
                method_name = col_methods.get(c, "First")

                # Internal merge methods
                if method_name in INTERNAL_AGG_FUNCS:
                    agg_row.append(INTERNAL_AGG_FUNCS[method_name](col_vals))
                    continue

                # Stats-based methods
                if method_name in STATS_METHODS:
                    numeric_vals = math_tool.remove_non_numeric(col_vals)
                    agg_row.append(stats_tool.descriptive_stat(numeric_vals, method_name, decimals))
                    continue

                raise ValueError(f"Unknown aggregation method: {method_name}")

            agg_by_key[key] = agg_row

        # Build output rows in original order
        out_rows = []
        emitted_keys = set()

        for row in data:
            key = row[dissolve_on]
            if key in agg_by_key:
                if key not in emitted_keys:
                    out_rows.append(agg_by_key[key])
                    emitted_keys.add(key)
            else:
                out_rows.append(row)

        # Build and return new table
        out_table = TableTools().new_table()
        out_table.set_data(out_rows)
        out_table.set_headers(headers)
        out_table.set_dtypes(dtypes)
        return out_table

    # ============================================================
    # TABLE-TO-TABLE OPERATIONS
    # ============================================================

    def join_table(self,to_join,primary_key,foreign_key,columns):
        """Join the data of this Table object with the data from another Table object ('to_join'), based on a primary key and foreign key. This join function is a left outer join, using a one-to-one relationship only (i.e., other join types and relationships are not supported).
        
        Parameters:
        
        to_join (object): the Table object to join to this Table object.
        
        primary_key (string): the column header of the column in this Table object that serves as the unique identifier.
        
        foreign_key (string): the column header of the column in the "to_join" object that corresponds with the primary key of this Table. Must be the same data type as the primary key.
        
        columns (list, string): a list of strings representing the column headers of the columns in the "to_join" data to be joined to this Table. May be specified as the string "all" to include all available columns."""
        if not self.column_exists(primary_key):
            raise ValueError(f"Primary key '{primary_key}' does not exist in this Table.")
        if not to_join.column_exists(foreign_key):
            raise ValueError(f"Foreign key '{foreign_key}' does not exist in to_join Table.")

        prim_key = self.get_column(primary_key)
        for_key = to_join.get_column(foreign_key)

        # Normalize columns argument
        if isinstance(columns, str):
            columns = [columns]
        if columns[0] == "all":
            columns = to_join.get_headers()

        # Validate requested columns
        for col_name in columns:
            if not to_join.column_exists(col_name):
                raise ValueError(f"Column '{col_name}' does not exist in to_join Table.")

        # Build lookup dict for foreign key 
        fk_lookup = {fk: i for i, fk in enumerate(for_key)}

        # For each requested column, build joined column
        for col_name in columns:
            col = to_join.get_column(col_name)
            out_col = []
            for pk in prim_key:
                if pk == "" or pk not in fk_lookup:
                    out_col.append("")
                else:
                    out_col.append(col[fk_lookup[pk]])
            self.append_column(out_col, col_name)
        return self

    def split_table(self,num_tables):
        """Split a Table into a specified number ('num_tables') of Table objects. Resulting Tables will have approximately equal numbers of rows, and will have the same headers and data types.
        
        Parameters:
        
        num_tables (integer): the number of Tables to split this Table into. Depending on the number of rows in the Table to be split, the actual number of split Tables may be one greater or lesser than the desired number."""
        lv = len(self._data)
        if lv == 0:
            return []

        if num_tables > lv:
            num_tables = lv

        sub_count = math.ceil(lv / num_tables)
        out_tables = []

        headers = self.get_headers()
        dtypes = self.get_dtypes()

        # Slice instead of popping
        for i in range(0, lv, sub_count):
            sub_list = self._data[i:i+sub_count]
            t = TableTools().new_table()
            t.set_data(sub_list)
            t.set_headers(headers)
            t.set_dtypes(dtypes)
            out_tables.append(t)

        return out_tables
    
    def append_table_rows(self,table):
        """Append the rows of another Table object ('Table') to the end of this Table. Tables must have the same number of columns. It is assumed column data types are the same between Tables.
        
        Parameters:
        
        table (object): the Table object to append rows from."""
        if self.is_empty():
            self.set_data(table.get_data())
            self.set_headers(table.get_headers())
            self.set_dtypes(table.get_dtypes())
        else:
            # Validate compatibility
            if self.num_columns() != table.num_columns():
                raise ValueError("Tables must have the same number of columns.")
            if self.get_dtypes() != table.get_dtypes():
                raise ValueError("Tables must have the same dtypes.")
            for r in range(table.num_rows()):
                self.append_row(table.get_row(r))
        return self
    
    def append_table_columns(self,table):
        """Append the columns of another Table object ('Table') to the end of this Table. Tables must have the same number of rows.
        
        Parameters:
        
        table (object): the Table object to append columns from."""
        if self.is_empty():
            self.set_data(table.get_data())
            self.set_headers(table.get_headers())
            self.set_dtypes(table.get_dtypes())
        else:
            if self.num_rows() != table.num_rows():
                raise ValueError("Tables must have the same number of rows.")
            else:
                other_headers = table.get_headers()
                for c in range(len(other_headers)):
                    col = table.get_column(c)
                    self.append_column(col,other_headers[c])
        return self

    def compare(self, table, summary=True):
        """Compare this Table with another Table object ('Table') and return a dictionary of structural differences. Returns a dictionary with keys: 'rows', 'columns', 'columns_only_in_self', 'columns_only_in_other', 'columns_in_both', 'column_type_differences', and 'column_order_matches'. Results are ordered self and other.

        Parameters:

        table (object): the Table object to compare against.
        
        summary (Boolean): flag to indicate if a summary of the comparison will be printed. If False, only returns a dictionary of results."""

        if not isinstance(table,_Table):
            raise ValueError("'table' must be a Table object.")

        rows_self = self.num_rows()
        rows_other = table.num_rows()
        cols_self = self.num_columns()
        cols_other = table.num_columns()

        headers_self = self.get_headers()
        headers_other = table.get_headers()

        dtypes_self = self.get_dtypes()
        dtypes_other = table.get_dtypes()

        set_self = set(headers_self)
        set_other = set(headers_other)

        only_self = sorted(list(set_self - set_other))
        only_other = sorted(list(set_other - set_self))
        both = sorted(list(set_self & set_other))

        type_diffs = {}
        for h in both:
            idx_self = headers_self.index(h)
            idx_other = headers_other.index(h)
            if dtypes_self[idx_self] != dtypes_other[idx_other]:
                type_diffs[h] = (dtypes_self[idx_self], dtypes_other[idx_other])

        order_matches = headers_self == headers_other

        result = {
            "rows": (rows_self, rows_other),
            "columns": (cols_self, cols_other),
            "columns_only_in_self": only_self,
            "columns_only_in_other": only_other,
            "columns_in_both": both,
            "column_type_differences": type_diffs,
            "column_order_matches": order_matches
        }

        if summary:
            print("Table Comparison Summary")
            print("------------------------")
            print("Rows:")
            print(f"  self:  {rows_self}")
            print(f"  other: {rows_other}")
            print("Columns:")
            print(f"  self:  {cols_self}")
            print(f"  other: {cols_other}")
            if only_self:
                print(f"Columns only in self: {only_self}")
            if only_other:
                print(f"Columns only in other: {only_other}")
            if both:
                print(f"Columns in both: {both}")
            if type_diffs:
                print("Column type differences:")
                for h, (t1, t2) in type_diffs.items():
                    print(f"  {h}: {t1} (self) vs {t2} (other)")
            print("Column order: " + ("matches" if order_matches else "differs"))
        return result

    def reshape(self,row=None,col=None):
        """Modify the number of rows and columns in the data by removing rows and columns beyond the specified row ('row') and column ('col') indices.
        
        Parameters:

        row (integer): the index of the last row to preserve. All rows beyond this index will be removed.

        col (integer): the index of the last column to preserve. All columns beyond this index will be removed."""
        if row is not None:
            if not isinstance(row, int) or row < 0 or row >= self.num_rows():
                raise ValueError("Row index is out of range.")
            # Apply row trimming
            self._data = self._data[:row + 1]

        # Validate col
        if col is not None:
            if not isinstance(col, int) or col < 0 or col >= self.num_columns():
                raise ValueError("Column index is out of range.")
            last_col = self.num_columns() - 1
            if col < last_col:
                self.remove_columns(col + 1, last_col)
        return self
            
    def transpose(self):
        """Transpose the data, flipping rows and columns and return a new Table."""
        data = self.get_data()  # safe copy

        if self.is_empty():
            raise ValueError("Cannot transpose an empty table.")

        # Transpose data matrix 
        transposed = list(map(list, zip(*data)))  # shape: num_columns x num_rows

        # Build output table
        _out = TableTools().new_table()
        _out.set_data(transposed)
        _out.autogenerate_headers()   # use your built-in header generator
        _out.detect_dtypes()          # re-infer dtypes for new orientation
        return _out

    def transpose_with_headers(self):
        """Transpose the data, flipping rows and columns, and include the original headers as the first column of the new table and return a new Table."""
        data = self.get_data()
        headers = self.get_headers()

        # Handle empty table 
        if self.is_empty():
            raise ValueError("Cannot transpose an empty table.")

        # Transpose data matrix
        transposed = list(map(list, zip(*data)))  # rows become original columns

        # Prepend original headers as first column values
        transposed_with_headers = [[headers[i]] + transposed[i] for i in range(self.num_columns())]

        # Build output table
        _out = TableTools().new_table()
        _out.set_data(transposed_with_headers)
        _out.autogenerate_headers()   # generate headers automatically
        _out.detect_dtypes()          # re-infer dtypes for new orientation
        return _out

    def pivot(self, index, columns, values):
        """Pivot the table from long to wide format and return a new Table.

        Parameters:

        index (integer, string): the index or header of the column to use as row identifiers.

        columns (integer, string): the index or header of the column whose unique values become new headers.

        values (integer, string): the index or header of the column whose values fill the new table."""
        # Normalize inputs
        for name, col in [("index", index), ("columns", columns), ("values", values)]:
            if isinstance(col, str):
                if not self.column_exists(col):
                    raise ValueError(f"Column '{col}' does not exist in Table.")
                locals()[name] = self._headers.index(col)
            elif isinstance(col, int):
                if col < 0 or col >= self.num_columns():
                    raise ValueError(f"Column index {col} is out of range.")
            else:
                raise TypeError(f"{name} must be a string or integer.")

        # Build dictionary keyed by index value
        pivot_dict = {}
        for row in self._data:
            key = row[index]
            col_name = row[columns]
            value = row[values]
            if key not in pivot_dict:
                pivot_dict[key] = {}
            pivot_dict[key][col_name] = value

        # Build output table
        new_headers = [self._headers[index]] + sorted({row[columns] for row in self._data})
        new_data = []
        for key, mapping in pivot_dict.items():
            new_data.append([key] + [mapping.get(h, "") for h in new_headers[1:]])

        _out = TableTools().new_table()
        _out.set_headers(new_headers)
        _out.set_data(new_data)
        _out.detect_dtypes()
        return _out

    def unpivot(self, id_vars, value_vars, var_name="Variable", value_name="Value"):
        """Unpivot the table from wide to long format and return a new Table.

        Parameters:

        id_vars (list): list of indices or headers to keep as identifier columns.

        value_vars (list): list of indices or headers to melt into rows.

        var_name (string): the header name for the new variable column.

        value_name (string): the header name for the new value column."""
        # Normalize id_vars
        id_indices = []
        for c in id_vars:
            if isinstance(c, str):
                if not self.column_exists(c):
                    raise ValueError(f"Column '{c}' does not exist in Table.")
                id_indices.append(self._headers.index(c))
            elif isinstance(c, int):
                if c < 0 or c >= self.num_columns():
                    raise ValueError(f"Column index {c} is out of range.")
                id_indices.append(c)
            else:
                raise TypeError("id_vars must contain only strings or integers.")

        # Normalize value_vars
        val_indices = []
        for c in value_vars:
            if isinstance(c, str):
                if not self.column_exists(c):
                    raise ValueError(f"Column '{c}' does not exist in Table.")
                val_indices.append(self._headers.index(c))
            elif isinstance(c, int):
                if c < 0 or c >= self.num_columns():
                    raise ValueError(f"Column index {c} is out of range.")
                val_indices.append(c)
            else:
                raise TypeError("value_vars must contain only strings or integers.")

        new_headers = [self._headers[i] for i in id_indices] + [var_name, value_name]
        new_data = []

        for row in self._data:
            id_part = [row[i] for i in id_indices]
            for vi in val_indices:
                val = row[vi] if row[vi] is not None else ""  # enforce empty string
                new_data.append(id_part + [self._headers[vi], val])

        _out = TableTools().new_table()
        _out.set_headers(new_headers)
        _out.set_data(new_data)
        _out.detect_dtypes()
        return _out
    
    def date_intersection(self,table,date_col1,date_col2):
        """Determine if the date column of another Table object ('table') overlaps with the date column of this Table and return both tables filtered by the overlapping range. It is assumed that the date ranges of each table are sorted in order. Returns a tuple of (filtered_table1, filtered_table2). If there is no overlap between the date ranges of the tables, a tuple of (None, None) is returned.
        
        Parameters:
        
        table (object): the Table object the operation will be performed on. The input Table will not be altered by the operation.

        date_col1 (integer, string): the index or header of the date column in this Table. Date column must be in a "yyyy-mm-dd" format. Any character may be used to separate date parts.

        date_col2 (integer, string): the index or header of the date column in the other Table. Date column must be in a "yyyy-mm-dd" format. Any character may be used to separate date parts."""
        if not isinstance(table,_Table):
            raise ValueError("'table' must be a Table object")
        if not self.column_exists(date_col1):
            raise ValueError(f"Column '{date_col1}' does not exist in this Table.")
        if not table.column_exists(date_col2):
            raise ValueError(f"Column '{date_col2}' does not exist in Table 'table'.")

        if isinstance(date_col1,str):
            date_col1 = self._headers.index(date_col1)
        if isinstance(date_col2,str):
            date_col2 = table._headers.index(date_col2)
        
        table1_start = self._data[0][date_col1]
        table1_end = self._data[-1][date_col1]

        table2_start = table._data[0][date_col2]
        table2_end = table._data[-1][date_col2]
        
        starts = [table1_start,table2_start]
        ends = [table1_end,table2_end]

        # Extract earliest/latest start/end.
        latest_start = max(starts)
        earliest_end = min(ends)

        # Now check for date range intersection.
        # If no date overlap, return None
        if (earliest_end < latest_start):
            return (None, None)
        
        else:
            # If there's perfect overlap (start and end dates are the same), we don't need to filter.
            if (table1_start == table2_start) and (table1_end == table2_end):
                return (self, table)
            
            # If there is partial overlap, only extract the overlap portions of the date range.
            elif (latest_start < earliest_end):
                out_table1 = self.copy()
                out_table2 = table.copy()
                if table1_start < latest_start:
                    out_table1 = out_table1.filter_date(date_col1,">=",latest_start)
                if table1_end > earliest_end:
                    out_table1 = out_table1.filter_date(date_col1,"<=",earliest_end)
                if table2_start < latest_start:
                    out_table2 = out_table2.filter_date(date_col2,">=",latest_start)
                if table2_end > earliest_end:
                    out_table2 = out_table2.filter_date(date_col2,"<=",earliest_end)
                return (out_table1, out_table2)

    # ============================================================
    # FILTERING
    # ============================================================

    def check_missing_invalid(self, missing_invalid = ""):
        """Check each column for missing or invalid values ('missing_invalid') and return a list of column indices that have missing or invalid values.
        
        Parameters:

        missing_invalid (integer, float, string, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value."""
        if not isinstance(missing_invalid, list):
            missing_invalid = [missing_invalid]

        output = []
        data = self.get_data()
        cols = self.num_columns()

        for col in range(cols):
            for row in data:
                if row[col] in missing_invalid:
                    output.append(col)
                    break  # no need to check further rows for this column

        return output

    def filter_numeric(self,column,condition,threshold):
        """Filter the Table data by removing rows that do not meet a specified threshold condition and return a new Table object. This filtering operation is performed with columns of numeric values.
        
        Parameters:
        
        column (integer, string): the index or header of the column to filter Table data by.
        
        condition (string): the condition placed upon each value in the specified column to determine inclusion in the filtered data. May be specified as "<" for less than, "<=" for less than equal to, ">" for greater than, ">=" for greater than equal to, "==" for equal to, or "!=" for not equal to.
        
        threshold (integer, float): the threshold to evaluate each value in the specified column against to determine inclusion in the filtered data."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in the table.")
        
        ops = {'<':operator.lt, '<=':operator.le, '>': operator.gt, '>=': operator.ge, '==': operator.eq, '!=': operator.ne}

        if condition not in ops:
            raise ValueError(
                f"Invalid condition '{condition}'. "
                f"Must be one of {list(ops.keys())}."
            )
        
        if not isinstance(threshold, (int, float)):
            raise TypeError("Threshold must be a numeric value.")

        if isinstance(column,str):
            column = self._headers.index(column)

        _out = TableTools().new_table()
        

        op = ops[condition]

        for row in self.get_data():
            val = row[column]
            if isinstance(val, (int, float)):
                if op(val, threshold):
                    _out.append_row(row[:])

        _out.set_headers(self.get_headers())
        _out.set_dtypes(self.get_dtypes())
        return _out
    
    def filter_text(self,column,condition,threshold):
        """Filter the Table data by removing rows that do not meet a specified threshold condition and return a new Table object. This filtering operation is performed with columns of string (text) values.
        
        Parameters:
        
        column (integer, string): the index or header of the column to filter Table data by.
        
        condition (string): the condition placed upon each value in the specified column to determine inclusion in the filtered data. May be specified as "==" for equal to, "!=" for not equal to, 'like' for contains, "not like" for does not contain, 'start' for starts with, 'not start' for does not start with, 'end' for ends with, or 'not end' for does not end with. Conditions are case sensitive.
        
        threshold (string): the threshold to evaluate each value in the specified column against to determine inclusion in the filtered data."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in the table.")
        
        valid_conditions = {
            "==",
            "!=",
            "like",
            "not like",
            "start",
            "not start",
            "end",
            "not end"
        }

        if condition not in valid_conditions:
            raise ValueError(
                f"Invalid condition '{condition}'. "
                f"Must be one of: {sorted(valid_conditions)}."
            )

        
        _out = TableTools().new_table()

        if isinstance(column,str):
            column = self._headers.index(column)

        for row in self.get_data():
            val = row[column]

            if not isinstance(val, str):
                continue  # non-string values fail the condition

            keep = False

            if condition == "==":
                keep = val == threshold
            elif condition == "!=":
                keep = val != threshold
            elif condition == "like":
                keep = threshold in val
            elif condition == "not like":
                keep = threshold not in val
            elif condition == "start":
                keep = val.startswith(threshold)
            elif condition == "not start":
                keep = not val.startswith(threshold)
            elif condition == "end":
                keep = val.endswith(threshold)
            elif condition == "not end":
                keep = not val.endswith(threshold)
            if keep:
                _out.append_row(row[:])

        _out.set_headers(self.get_headers())
        _out.set_dtypes(self.get_dtypes())
        return _out

    def filter_date(self,column,condition,threshold):
        """Filter the Table data by removing rows that do not meet a specified threshold condition and return a new Table object. This filtering operation is performed with columns of date values, which must be in "yyyy-mm-dd" format.
        
        Parameters:
        
        column (integer, string): the index or header of the column to filter Table data by. Date values in this column must be in "yyyy-mm-dd" format
        
        condition (string): the condition placed upon each value in the specified column to determine inclusion in the filtered data. May be specified as "<" for less than, "<=" for less than equal to, ">" for greater than, ">=" for greater than equal to, "==" for equal to, or "!=" for not equal to.
        
        threshold (string): the threshold to evaluate each value in the specified column against to determine inclusion in the filtered data. Must be in "yyyy-mm-dd" format."""
        # Validate column existence
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in the table.")

        # Supported operators
        ops = {
            '<': operator.lt,
            '<=': operator.le,
            '>': operator.gt,
            '>=': operator.ge,
            '==': operator.eq,
            '!=': operator.ne
        }

        if condition not in ops:
            raise ValueError(
                f"Invalid condition '{condition}'. "
                f"Must be one of {list(ops.keys())}."
            )
        
        if not isinstance(threshold, str) or len(threshold) != 10:
            raise ValueError("Threshold must be a date string in 'YYYY-MM-DD' format.")

        _out = TableTools().new_table()

        if isinstance(column,str):
            column = self._headers.index(column)

        op = ops[condition]

        for row in self.get_data():
            val = row[column]

            # Only evaluate ISO date strings
            if isinstance(val, str) and len(val) == 10:
                if op(val, threshold):
                    _out.append_row(row[:])

        _out.set_headers(self.get_headers())
        _out.set_dtypes(self.get_dtypes())
        return _out

    def filter_boolean(self, column, condition, target):
        """Filter the Table data by removing rows that do not meet a specified Boolean condition and return a new Table object.

        Parameters:

        column (integer, string): the index or header of the column to filter by.

        condition (string): the condition used to evaluate each Boolean value. Must be "==" for equality or "!=" for inequality.

        target (bool): the Boolean value to compare against."""

        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in the table.")

        # Only equality and inequality make sense for booleans
        ops = {"==": operator.eq, "!=": operator.ne}

        if condition not in ops:
            raise ValueError(
                f"Invalid condition '{condition}'. "
                f"Must be one of {list(ops.keys())}."
            )

        if not isinstance(target, bool):
            raise TypeError("Target must be a Boolean value (True or False).")

        if isinstance(column, str):
            column = self._headers.index(column)

        _out = TableTools().new_table()
        op = ops[condition]

        for row in self.get_data():
            val = row[column]
            if isinstance(val, bool):
                if op(val, target):
                    _out.append_row(row[:])

        _out.set_headers(self.get_headers())
        _out.set_dtypes(self.get_dtypes())
        return _out

    def filter_type(self, column, condition, dtype):
        """Filter the Table data by removing rows that do not meet a specified data type condition and return a new Table object. This filtering operation is performed with columns of specific data types.
        
        Parameters:
        
        column (integer, string): the index or header of the column to filter Table data by.
        
        condition (string): the condition placed upon each value in the specified column to determine inclusion in the filtered data. May be specified as "is" to include values that match the data type, or "is not" to include values that do not match.
        
        dtype (string): the data type to evaluate each value in the specified column against to determine inclusion in the filtered data. Possible types include "int", "float", "string", "bool", and "nan"."""
        # Validate column existence
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in the table.")

        # Resolve column index
        if isinstance(column, str):
            column = self._headers.index(column)

        # Validate condition
        if condition not in ("is", "is not"):
            raise ValueError("Condition must be 'is' or 'is not'.")

        # Validate dtype
        valid_types = {"int", "float", "string", "bool", "nan"}
        if dtype not in valid_types:
            raise ValueError(f"Invalid dtype '{dtype}'. Must be one of {valid_types}.")

        def matches(val):
            if dtype == "int":
                # Exclude bools since bool is a subclass of int
                return isinstance(val, int) and not isinstance(val, bool)
            elif dtype == "float":
                return isinstance(val, float) and not math.isnan(val)
            elif dtype == "string":
                return isinstance(val, str)
            elif dtype == "bool":
                return isinstance(val, bool)
            elif dtype == "nan":
                return isinstance(val, float) and math.isnan(val)
            return False

        _out = TableTools().new_table()

        for row in self.get_data():
            val = row[column]
            is_match = matches(val)
            keep = (is_match if condition == "is" else not is_match)
            if keep:
                _out.append_row(row[:])

        # Set headers and dtypes AFTER rows are added
        _out.set_headers(self.get_headers())
        _out.set_dtypes(self.get_dtypes())

        return _out

    def filter_missing_invalid(self,to_remove = ""):
        """Filter the Table data by removing rows that contain missing or invalid values ('to_remove') and return a new Table.
        
        Parameters:

        to_remove (integer, float, string, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value."""
        if not isinstance(to_remove, list):
            to_remove = [to_remove]

        _out = TableTools().new_table()

        for row in self.get_data():
            # Skip rows containing any invalid value
            if any(val in to_remove for val in row):
                continue
            _out.append_row(row[:])

        # Preserve headers and dtypes
        _out.set_headers(self.get_headers())
        _out.set_dtypes(self.get_dtypes())

        return _out
  
    # ============================================================
    # CALCULATIONS AND EXPRESSIONS
    # ============================================================

    def generate_id_col(self,header = "ID", start = 1, step = 1, index = 0):
        """Generate an ID column of unique integer values.
        
        Parameters:

        header (string): the header of the ID column.

        start (integer): the first value of the ID column.

        step (integer): the difference between subsequent values. May be positive (for addition) or negative (for subtraction).

        index (integer): the index where the column will be inserted."""
        n = self.num_rows()

        if n == 0:
            raise ValueError("Cannot generate an ID column for an empty table.")

        # Convert header to string
        header = str(header)

        # Validate index
        if not isinstance(index, int):
            raise TypeError("Index must be an integer.")

        if index < 0 or index > self.num_columns():
            raise ValueError(f"Index {index} is out of range for insertion.")

        # Generate the sequence
        output = [start + i * step for i in range(n)]

        # Insert the column
        self.insert_column(output, index, header)
        return self
   
    def column_calculator(self, expression, column=None, missing_invalid="", in_place=True, new_header=None, decimals=2):
        """Evaluate a specified expression ('expression') for each row in the Table object using column references, row-level functions such as mathematical and conditional operators, column-level scalar, and column transformers. The results of the calculator can overwrite a target column or be appended as a new column. The helper function calculator_help() may be called to explain calculator syntax and available functions.

        Parameters:
        
        expression (string): the expression to evaluate. 
        
        column (integer, string): the target column the calculator results will be written to if the in_place flag is True.

        missing_invalid (string, list): a string or list of values specifying missing or invalid values in the Table. If a missing or invalid value is present in a row, the column calculator will skip that row.
        
        in_place (Boolean): flag to indicate if the results of the column calculator will overwrite the target column. Calculator results will be appended as a new column if flag is False.
        
        new_header (string): if in_place flag is False, the header of the new column to be appended.
        
        decimals (integer): the number of decimal places float values will be rounded to."""
        if in_place and column is None:
            raise ValueError("When in_place=True, you must specify a target column to overwrite.")

        if in_place:
            if not self.column_exists(column):
                raise ValueError(f"Column '{column}' does not exist in Table.")

            if isinstance(column, str):
                col_index = self._headers.index(column)
                col_name = column
            else:
                col_index = column
                col_name = self._headers[column]
        else:
            # Appending a new column — no target column needed
            col_index = None
            col_name = None

        # Normalize missing/invalid list
        if not isinstance(missing_invalid, list):
            missing_invalid = [missing_invalid]

        # Build safe environment for eval
        safe_env = {"__builtins__": {}}
        safe_env.update(self._SAFE_ENV)

        # Regex patterns for parsing
        import re
        import ast

        # Matches {Header}
        header_pattern = re.compile(r"\{([^}]+)\}")

        # Matches col@#
        index_pattern = re.compile(r"col@(\d+)")

        # Matches AGGREGATE_FUNCTION(args) — ALL CAPS
        agg_pattern = re.compile(r"([A-Z][A-Z0-9_]*)\s*\((.*?)\)")

        # Matches lowercase function calls: func(...)
        # This is for COLUMN_FUNCTIONS (not row-level builtins)
        colfunc_pattern = re.compile(r"([a-z_][a-z0-9_]*)\s*\((.*?)\)")

        # Extract referenced headers and indices (initial scan)
        referenced_headers = header_pattern.findall(expression)
        referenced_indices = [int(i) for i in index_pattern.findall(expression)]

        # Validate referenced headers
        for h in referenced_headers:
            if h not in self._headers:
                raise ValueError(f"Expression references unknown column header: '{h}'")

        # Validate referenced indices
        for idx in referenced_indices:
            if not self.column_exists(idx):
                raise ValueError(f"Expression references invalid column index: col@{idx}")

        # Detect column-level transformer functions
    
        # These are lowercase functions defined in self._COLUMN_FUNCTIONS.
        # They operate on entire columns and return full lists.
        # We detect them BEFORE aggregates and BEFORE row-level eval.
        # If found, we evaluate them recursively and replace them with
        # placeholder variables in the expression.

        # This dictionary will map placeholder names to full column lists
        columnfunc_results = {}

        # We will fill this in Chunk 2
        # (placeholder: actual logic added in next chunk)

        # COLUMN-FUNCTION DETECTION & EVALUATION
        # Column-functions: are lowercase, operate on entire columns, return full lists, must be evaluated BEFORE aggregates and row-level eval
        # We recursively detect and evaluate any column-function calls
        # in the expression, replacing them with placeholder variables.
        def resolve_argument(arg_text):
            """
            Resolve a single argument for a column-function or scalar:
            • {Header}  : full column list
            • col@#     : full column list
            • literal   : parsed via ast.literal_eval
            """
            # Column reference by header
            header_match = header_pattern.fullmatch(arg_text)
            if header_match:
                header = header_match.group(1)
                if header not in self._headers:
                    raise ValueError(f"Unknown column header '{header}' in function call.")
                return self.get_column(header)

            # Column reference by index
            index_match = index_pattern.fullmatch(arg_text)
            if index_match:
                idx = int(index_match.group(1))
                if not self.column_exists(idx):
                    raise ValueError(f"Invalid column index col@{idx} in function call.")
                return self.get_column(idx)

            # Literal value (string, number, bool, None)
            try:
                return ast.literal_eval(arg_text)
            except Exception:
                raise ValueError(f"Invalid argument '{arg_text}'. Must be a column reference or literal.")

        # Recursive evaluation of column-functions
        def evaluate_column_functions(expr):
            """
            Detects and evaluates column-functions in the expression.
            Returns:
                new_expr  : expression with placeholders inserted
                results   : dict mapping placeholder to list
            """
            results = {}
            changed = True
            placeholder_counter = 0

            while changed:
                changed = False
                matches = colfunc_pattern.findall(expr)

                for func_name, arg_string in matches:
                    # Skip if not a registered column-function
                    if func_name not in self._COLUMN_FUNCTIONS:
                        continue

                    changed = True

                    # Split argument list
                    raw_args = [a.strip() for a in arg_string.split(",") if a.strip()]
                    if len(raw_args) == 0:
                        raise ValueError(f"Column-function '{func_name}' requires at least one argument.")

                    # Resolve arguments (columns or literals)
                    resolved_args = [resolve_argument(a) for a in raw_args]

                    # Evaluate the column-function
                    try:
                        import inspect

                        func = self._COLUMN_FUNCTIONS[func_name]["func"]
                        sig = inspect.signature(func)
                        params = sig.parameters

                        # Does the function explicitly accept a "decimals" parameter?
                        accepts_decimals = "decimals" in params

                        # Does the function accept **kwargs?
                        accepts_kwargs = any(
                            p.kind == inspect.Parameter.VAR_KEYWORD
                            for p in params.values()
                        )

                        if accepts_decimals or accepts_kwargs:
                            result_list = func(*resolved_args, decimals=decimals)
                        else:
                            result_list = func(*resolved_args)

                    except Exception as e:
                        raise ValueError(f"Error computing column-function '{func_name}': {e}")

                    # Create placeholder
                    placeholder = f"_COLFUNC_{func_name}_{placeholder_counter}"
                    placeholder_counter += 1

                    # Store result list
                    results[placeholder] = result_list

                    # Replace the function call in the expression
                    full_call_text = f"{func_name}({arg_string})"
                    expr = expr.replace(full_call_text, placeholder)

            return expr, results

        # Apply column-function evaluation to the expression
        rewritten_expr, columnfunc_results = evaluate_column_functions(expression)

        # Add column-function results to safe_env
        
        # DETECT AND COMPUTE AGGREGATE FUNCTIONS (ALL CAPS)
        # Scalar: operate on entire columns, return a single scalar, evaluated ONCE, replaced with placeholder variables

        # Find all aggregate calls in the rewritten expression
        aggregate_calls = agg_pattern.findall(rewritten_expr)

        # This will store mappings from original text to safe variable name
        agg_replacements = {}

        # Counter to ensure unique safe variable names
        agg_counter = 0

        for func_name, arg_string in aggregate_calls:

            # Skip if function is not in our registry
            if func_name not in self._SCALAR_FUNCTIONS:
                raise ValueError(f"Unknown scalar function '{func_name}' in expression.")

            # Parse argument list inside parentheses
            raw_args = [arg.strip() for arg in arg_string.split(",") if arg.strip()]

            if len(raw_args) == 0:
                raise ValueError(f"Scalar function '{func_name}' requires at least one argument.")

            # Resolve each argument to a full column list
            # (Aggregates still only accept column references)
            col_lists = []

            for arg in raw_args:

                # Case 1: {Header}
                header_match = header_pattern.fullmatch(arg)
                if header_match:
                    header = header_match.group(1)
                    if header not in self._headers:
                        raise ValueError(f"Unknown column header '{header}' in scalar call.")
                    col_lists.append(self.get_column(header))
                    continue

                # Case 2: col@#
                index_match = index_pattern.fullmatch(arg)
                if index_match:
                    idx = int(index_match.group(1))
                    if not self.column_exists(idx):
                        raise ValueError(f"Invalid column index col@{idx} in scalar call.")
                    col_lists.append(self.get_column(idx))
                    continue

                # If we reach here, the argument is invalid
                raise ValueError(
                    f"Invalid argument '{arg}' in scalar function '{func_name}'. "
                    "Scalar arguments must be {Header} or col@#."
                )

            # Compute the aggregate value using the registry
            try:
                agg_value = self._SCALAR_FUNCTIONS[func_name]["func"](*col_lists, decimals)
            except Exception as e:
                raise ValueError(
                    f"Error computing scalar '{func_name}': {e}"
                )

            # Create a safe variable name for this aggregate
            safe_agg_name = f"_AGG_{func_name}_{agg_counter}"
            agg_counter += 1

            safe_env[safe_agg_name] = agg_value

            # Store replacement mapping
            full_call_text = f"{func_name}({arg_string})"
            agg_replacements[full_call_text] = safe_agg_name

        # Rewrite expression to replace aggregate calls with safe names
        for original, safe_name in agg_replacements.items():
            rewritten_expr = rewritten_expr.replace(original, safe_name)

        # Rewrite {Header} and col@# references for row-level eval
        # At this point, rewritten_expr has: column-function placeholders and aggregate placeholders
        # We now map: {Header} to safe per-row variable names, col@# to safe per-row variable names

        # Extract headers and indices again from the rewritten expression
        referenced_headers = header_pattern.findall(rewritten_expr)
        referenced_indices = [int(i) for i in index_pattern.findall(rewritten_expr)]

        # Map headers to safe variable names
        header_var_map = {}
        for h in referenced_headers:
            if h not in header_var_map:
                safe_name = f"_col_header_{self._headers.index(h)}"
                header_var_map[h] = safe_name
                rewritten_expr = rewritten_expr.replace(f"{{{h}}}", safe_name)

        # Map indices to safe variable names
        index_var_map = {}
        for idx in referenced_indices:
            if idx not in index_var_map:
                safe_name = f"_col_index_{idx}"
                index_var_map[idx] = safe_name
                rewritten_expr = rewritten_expr.replace(f"col@{idx}", safe_name)

        new_vals = []
        success = True

        # Row-by-row evaluation
        for r in range(self.num_rows()):
            row = self.get_row(r)

            # Build row variable dictionary
            row_vars = {}

            # Add header-based variables
            for h, safe_name in header_var_map.items():
                val = row[self._headers.index(h)]
                row_vars[safe_name] = val

            # Add index-based variables
            for idx, safe_name in index_var_map.items():
                val = row[idx]
                row_vars[safe_name] = val

            # Add column-function per-row values
            for placeholder, col_list in columnfunc_results.items():
                # col_list is the full column produced by the column-function
                # We take the value for the current row r
                row_vars[placeholder] = col_list[r]

            # Check for missing/invalid values (now includes column-function values)
            if any(v in missing_invalid for v in row_vars.values()):
                # Preserve original value if skipping
                new_vals.append(self._data[r][col_index])
                continue

            # Evaluate safely
            try:
                # Combine safe_env (aggregates + math + column-functions)
                # with row_vars (per-row values)
                result = eval(rewritten_expr, safe_env, row_vars)

                # Apply rounding if result is float
                if isinstance(result, float):
                    result = round(result, decimals)

                new_vals.append(result)

            except Exception as e:
                success = False
                print(f"Error in Column Calculator in row {r+1}: {e}")
                break
        
        # Write results back to the table
        if success:
            if in_place:
                # Overwrite the existing column
                self.replace_column(col_index, new_vals, col_name)
            else:
                # Append as a new column
                if new_header == "":
                    raise ValueError("Please provide a header for the new column.")
                self.append_column(new_vals, new_header)

        return self

    def calculator_help(self, terminal=False):
        """Display help for the column calculator. Generates an HTML document by default.

        Parameters:

        terminal (Boolean): flag to indicate if the help will be printed to a terminal or written as a file."""
        random.seed(12345)  # stable examples

        # Filter out __builtins__ from row-level functions
        safe_env_keys = [k for k in self._SAFE_ENV.keys() if k != "__builtins__"]

        # Helper: choose up to 3 example functions from a registry
        def pick_examples(keys):
            if not keys:
                return []
            return random.sample(keys, k=min(3, len(keys)))

        # Pick examples dynamically
        row_examples = pick_examples(safe_env_keys)
        agg_examples = pick_examples(list(self._SCALAR_FUNCTIONS.keys()))
        col_examples = pick_examples(list(self._COLUMN_FUNCTIONS.keys()))

        # TERMINAL MODE
        if terminal:
            print("TABLETOOLS COLUMN CALCULATOR")
            print("The TableTools Column Calculator lets you create or overwrite columns "
                "by building expressions that operate on your data using Python‑style syntax, "
                "row‑level functions, column‑level scalars, and full‑column transformer functions.\n")
            print("This help page lists all available functions, shows how to reference columns, "
                "and provides examples drawn directly from the functions currently loaded in your "
                "TableTools installation.\n")

            print("COLUMN REFERENCES")
            print("  {Header}     Reference a column by header")
            print("  col@#        Reference a column by index\n")

            print("ROW-LEVEL FUNCTIONS (evaluated per row)")
            for name in sorted(safe_env_keys):
                print(f"    • {name}")
            if row_examples:
                print("\nExamples:")
                for fn in row_examples:
                    print(f"  {fn}({{Flow}})")
            print()

            print("COLUMN-LEVEL SCALAR FUNCTIONS (ALL CAPS)")
            for key in sorted(self._SCALAR_FUNCTIONS.keys()):
                meta = self._SCALAR_FUNCTIONS[key]
                print(f"    • {key:<8}  {meta.get('name','')}")
            if agg_examples:
                print("\nExamples:")
                for fn in agg_examples:
                    print(f"  {fn}({{Flow}})")
            print()

            print("COLUMN-LEVEL TRANSFORMER FUNCTIONS (lowercase)")
            for key in sorted(self._COLUMN_FUNCTIONS.keys()):
                meta = self._COLUMN_FUNCTIONS[key]
                print(f"    • {key:<12}  {meta.get('name','')}")
            if col_examples:
                print("\nExamples:")
                for fn in col_examples:
                    print(f"  {fn}({{Flow}})")
            print()

            print("FULL EXPRESSION EXAMPLES")
            print("  ({Flow} - MEAN({Flow})) / STDEV({Flow})")
            print("  math.sqrt({X}**2 + {Y}**2)")
            print("  round({Flow} * 1.2, 3)")
            print()
            return

        # HTML MODE
        import webbrowser as wb

        package_dir = os.path.dirname(os.path.abspath(__file__))
        docs_dir = os.path.join(package_dir, "Docs")
        html_path = os.path.join(docs_dir, "calculator_help.html")

        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Build dynamic lists
        row_funcs_html = "\n".join(
            "<li><code>" + name + "</code></li>"
            for name in sorted(safe_env_keys)
        )
        agg_funcs_html = "\n".join(
            "<li><code>" + key + "</code> — " +
            self._SCALAR_FUNCTIONS[key].get("name", "") +
            "</li>"
            for key in sorted(self._SCALAR_FUNCTIONS.keys())
        )
        col_funcs_html = "\n".join(
            "<li><code>" + key + "</code> — " +
            self._COLUMN_FUNCTIONS[key].get("name", "") +
            "</li>"
            for key in sorted(self._COLUMN_FUNCTIONS.keys())
        )

        # Example blocks — always use {Flow}
        def example_block(funcs):
            if not funcs:
                return ""
            return "\n".join("  " + fn + "({Flow})" for fn in funcs)

        row_example_html = example_block(row_examples)
        agg_example_html = example_block(agg_examples)
        col_example_html = example_block(col_examples)

        html = """
        <!DOCTYPE html>
        <html>
        <head>
        <meta charset="UTF-8">
        <title>TableTools Column Calculator Help</title>
        <link rel="stylesheet" href="styles.css">
        </head>

        <body>
        <div class="container">

            <div class="sidebar">
                <div class="side-head">Column Calculator Help</div>
                <a href="#preamble" class="func">Overview</a>
                <a href="#colrefs" class="func">Column references</a>
                <a href="#rowfuncs" class="func">Row-level functions</a>
                <a href="#scalfuncs" class="func">Scalar functions</a>
                <a href="#colfuncs" class="func">Column-level functions</a>
                <a href="#examples" class="func">Full examples</a>
            </div>

            <div class="main-content">
                <h1>TableTools Column Calculator Help</h1>
                <p class="updated-timestamp">Generated: """ + timestamp + """</p>

                <h2 id="preamble">Overview</h2>
                <div class="function-doc">
                    <p>The TableTools Column Calculator lets you create or overwrite columns by building expressions that operate on your data using Python‑style syntax, row‑level functions, column‑level scalars, and full‑column transformer functions.</p>
                    <p>This help page lists all available functions, shows how to reference columns, and provides examples drawn directly from the functions currently loaded in your TableTools installation.</p>
                </div>

                <h2 id="colrefs">Column References</h2>
                <div class="function-doc">
                    <p><code>{Header}</code> — reference a column by header</p>
                    <p><code>col@#</code> — reference a column by index</p>
                    <pre>
        Examples:
        {Flow}
        col@2
        ({X}**2 + {Y}**2)**0.5
                    </pre>
                </div>

                <h2 id="rowfuncs">Row-Level Functions</h2>
                <div class="function-doc">
                    <ul>
        """ + row_funcs_html + """
                    </ul>
        """ + (("<pre>\n" + row_example_html + "\n</pre>") if row_example_html else "") + """
                </div>

                <h2 id="scalfuncs">Scalar Functions (ALL CAPS)</h2>
                <div class="function-doc">
                    <ul>
        """ + agg_funcs_html + """
                    </ul>
        """ + (("<pre>\n" + agg_example_html + "\n</pre>") if agg_example_html else "") + """
                </div>

                <h2 id="colfuncs">Column-Level Transformer Functions (lowercase)</h2>
                <div class="function-doc">
                    <ul>
        """ + col_funcs_html + """
                    </ul>
        """ + (("<pre>\n" + col_example_html + "\n</pre>") if col_example_html else "") + """
                </div>

                <h2 id="examples">Full Expression Examples</h2>
                <div class="function-doc">
                    <pre>
        ({Flow} - MEAN({Flow})) / STDEV({Flow})
        math.sqrt({X}**2 + {Y}**2)
        round({Flow} * 1.2, 3)
                    </pre>
                </div>

            </div>
        </div>
        </body>
        </html>"""

        with open(html_path, "w", encoding="utf-8") as f:
            f.write(html)

        wb.open_new_tab(html_path)

    def apply_column(self, column, func, args = None, kwargs = None, new_header=None, return_result=False):
        """Apply a user-defined function element-wise to a column ('column'). The function ('func') is applied to each value in the specified column. The result can be returned directly as a list, used to overwrite the existing column, or appended as a new column.

        Parameters:

        column (integer, string): the index or header of the column the operation will be performed on.

        func (callable): the function that accepts a single value and returns a transformed value. The function is applied element-wise.

        args (tuple): the tuple of additional positional arguments to pass to the user-defined function. No lambda expressions are required.

        kwargs (dict): the dictionary of additional keyword arguments to pass to the user-defined function. No lambda expressions are required.

        new_header (string, optional): if provided and return_result is False, the transformed values will be appended as a new column with this header.

        return_result (bool): if True, the transformed list is returned and the table is not modified."""
        if not callable(func):
            raise TypeError("func must be a callable function.")
        if args is not None and not isinstance(args, tuple):
            raise TypeError("args must be a tuple of positional arguments.")
        if kwargs is not None and not isinstance(kwargs, dict):
            raise TypeError("kwargs must be a dictionary of keyword arguments.")
        if new_header is not None:
            if not isinstance(new_header, str) or new_header == "":
                raise TypeError("new_header must be a non-empty string.")
            if new_header in self._headers:
                raise ValueError(f"Column header '{new_header}' already exists in the table.")
        if not isinstance(return_result, bool):
            raise TypeError("return_result must be a boolean.")
        # Get the original column values (validated internally)
        original = self.get_column(column)

        # normalize arguments, if any
        if args is None:
            args = ()
        if kwargs is None:
            kwargs = {}

        # apply function
        results = [func(v, *args, **kwargs) for v in original]

        # If returning results, do not modify the table
        if return_result:
            if new_header is not None:
                raise ValueError(f"{new_header} cannot be appended as a new column when when return_result=True")
            return results

        # Overwrite the existing column
        if new_header is None:
            self.replace_column(column, results, self._headers[self._resolve_column_index(column)])
            return self

        # Append as a new column
        self.append_column(results, new_header)
        return self

    def map_column(self, column, func, args=None, kwargs=None):
        """Apply a user-defined function ('func') to an entire column ('column') and return the result.

        The function ('func') receives the full list of column values as its first argument. The function may return any type of object, including a list of any length, a scalar, or a complex structure. The Table is not modified.

        Parameters:

        column (integer, string): the index or header of the column to be passed to the function.

        func (callable): the function that accepts a list of values and returns any result.

        args (tuple): additional positional arguments to pass to the user-defined function. No lambda expressions are required.

        kwargs (dict, optional): additional keyword arguments to pass to the user-defined function. No lambda expressions are required."""
        if not callable(func):
            raise TypeError("func must be a callable function.")
        if args is not None and not isinstance(args, tuple):
            raise TypeError("args must be a tuple of positional arguments.")
        if kwargs is not None and not isinstance(kwargs, dict):
            raise TypeError("kwargs must be a dictionary of keyword arguments.")
        
        # Normalize args/kwargs
        if args is None:
            args = ()
        if kwargs is None:
            kwargs = {}

        # Get the full column list
        col = self.get_column(column)

        # Apply the function to the entire list
        return func(col, *args, **kwargs)

    def apply_row(self, row, func, args=None, kwargs=None, append_row=False, return_result=False):
        """Apply a user-defined function element-wise to a row ('row'). The function ('func') is applied to each value in the specified row. The result can be returned directly as a list, used to overwrite the existing row, or appended as a new row.

        Parameters:

        row (integer): the index of the row the operation will be performed on.

        func (callable): the function that accepts a single value and returns a transformed value. The function is applied element-wise.

        args (tuple): the tuple of additional positional arguments to pass to the user-defined function. No lambda expressions are required.

        kwargs (dict): the dictionary of additional keyword arguments to pass to the user-defined function. No lambda expressions are required.

        append_row (bool): if True and return_result is False, the transformed row will be appended as a new row at the end of the table. If False and return_result is False, the transformed row will overwrite the existing row at 'row_index'."""
        # Validate func
        if not callable(func):
            raise TypeError("func must be a callable function.")
        # Validate args
        if args is not None and not isinstance(args, tuple):
            raise TypeError("args must be a tuple of positional arguments.")
        # Validate kwargs
        if kwargs is not None and not isinstance(kwargs, dict):
            raise TypeError("kwargs must be a dictionary of keyword arguments.")
        # Validate append_row and return_result
        if not isinstance(append_row, bool):
            raise TypeError("append_row must be a boolean.")
        if not isinstance(return_result, bool):
            raise TypeError("return_result must be a boolean.")
        if append_row and return_result:
                raise ValueError("row cannot be appended as a new row when when return_result=True")

        # Get the original row values (validated internally)
        original = self.get_row(row)

        # Normalize arguments, if any
        if args is None:
            args = ()
        if kwargs is None:
            kwargs = {}

        # Apply function element-wise
        results = [func(v, *args, **kwargs) for v in original]

        # If returning results, do not modify the table
        if return_result:
            return results

        # Overwrite the existing row or append as a new row
        if append_row:
            self.append_row(results)
        else:
            self.replace_row(row, results)

        return self

    def map_row(self, row, func, args=None, kwargs=None):
        """Apply a user-defined function ('func') to an entire row ('row') and return the result. The function ('func') receives the full list of row values as its first argument. The function may return any type of object, including a list of any length, a scalar, or a complex structure. The Table is not modified.

        Parameters:

        row (integer): the index of the row to be passed to the function.

        func (callable): the function that accepts a list of values and returns any result.

        args (tuple): additional positional arguments to pass to the user-defined function. No lambda expressions are required.

        kwargs (dict, optional): additional keyword arguments to pass to the user-defined function. No lambda expressions are required."""
        # Validate func
        if not callable(func):
            raise TypeError("func must be a callable function.")

        # Validate args
        if args is not None and not isinstance(args, tuple):
            raise TypeError("args must be a tuple of positional arguments.")

        # Validate kwargs
        if kwargs is not None and not isinstance(kwargs, dict):
            raise TypeError("kwargs must be a dictionary of keyword arguments.")

        # Normalize args/kwargs
        if args is None:
            args = ()
        if kwargs is None:
            kwargs = {}

        # Get the full row list
        row = self.get_row(row)

        # Apply the function to the entire row
        return func(row, *args, **kwargs)

    def round_column(self,column,decimals = 2):
        """Round the values in a specified column ('column') to a specified number of decimal places ('decimals'). Only numeric values (integer or float) in the specified column will be rounded.
        
        Parameters:

        column (integer, string): the index or header of the column the operation will be performed on.
        
        decimals (integer): the number of decimal places float values will be rounded to."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")
        else:
            if isinstance(column,str):
                column = self._headers.index(column)

            for i, row in enumerate(self._data):
                value = row[column]
                if isinstance(value, (int, float)):
                    self._data[i][column] = round(value,decimals)
        return self

    def column_stat(self,column,stat,decimals = 2):
        """Calculate and return the specified descriptive statistic ('stat') of the specified column ('column'). Only numeric values (integer or float) in the specified column will be included in the calculation.

        Parameters:

        column (integer, string): the index or header of the column the operation will be performed on.

        stat (string): the specific statistic to return. Possible statistics are "Count", "Unique Values", "First", "Last", "Sum", "Max", "Min", "Range", "Mean", "Median", "Mode", "IQ1", "IQ3", "IQR", "Standard Deviation", "Standard Error", "Variance", "Coefficient of Variation", "Skewness", and "Kurtosis"
        
        decimals (integer): the number of decimal places float values will be rounded to."""
        
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")
        
        valid_stats = {
        "Count", "Unique Values", "First", "Last", "Sum", "Max", "Min", "Range", "Mean",
        "Median", "Mode", "IQ1", "IQ3", "IQR", "Standard Deviation",
        "Standard Error", "Variance", "Coefficient of Variation",
        "Skewness", "Kurtosis"}

        if stat not in valid_stats:
            raise ValueError(f"Unknown statistic '{stat}'.")
        
        c = self.get_column(column)
        c = _Math().remove_non_numeric(c)
        if not c:
            return None
        return _Stats().descriptive_stat(c,stat,decimals)
    
    def column_summary(self,column,terminal = False, decimals = 2):
        """Calculate all descriptive statistics of the specified column ('column') and return a dictionary of all statistics. Only numeric values (integer or float) in the specified column will be included in the calculation.

        Parameters:

        column (integer, string): the index or header of the column the operation will be performed on.

        terminal (Boolean): flag to indicate if the result will be printed to a terminal or just returned.

        decimals (integer): the number of decimal places float values will be rounded to."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")
        if not isinstance(column,str):
            column = self._headers[column]
        
        c = self.get_column(column)
        c = _Math().remove_non_numeric(c)
        if not c:
            return None
        return _Stats().summary_stats(c, terminal, decimals)

    def column_stats_list(self,stat,decimals = 2):
        """Calculate the specified descriptive statistic ('stat') of each column in the Table and return a list of statistics. Only numeric values (integer or float) in each column will be included in the calculation.

        Parameters:

        stat (string): the specific statistic to calculate for each column. Possible statistics are "Count", "Unique Values", "First", "Last", "Sum", "Max", "Min", "Range", "Mean", "Median", "Mode", "IQ1", "IQ3", "IQR", "Standard Deviation", "Standard Error", "Variance", "Coefficient of Variation", "Skewness", and "Kurtosis"
        
        decimals (integer): the number of decimal places float values will be rounded to."""
        valid_stats = {
        "Count", "Unique Values", "First", "Last", "Sum", "Max", "Min", "Range", "Mean",
        "Median", "Mode", "IQ1", "IQ3", "IQR", "Standard Deviation",
        "Standard Error", "Variance", "Coefficient of Variation",
        "Skewness", "Kurtosis"}

        if stat not in valid_stats:
            raise ValueError(f"Unknown statistic '{stat}'.")
        
        output = []
        stats = _Stats()
        maths = _Math()
        for i in range(self.num_columns()):
            c = self.get_column(i)
            c = maths.remove_non_numeric(c)
            if not c:
                output.append(None)
            else:
                output.append(stats.descriptive_stat(c,stat,decimals))
        return output

    def row_stat(self,row,stat,decimals = 2):
        """Calculate and return the specified descriptive statistic ('stat') of the specified row ('row'). Only numeric values (integer or float) in the specified row will be included in the calculation.

        Parameters:

        row (integer): the index of the row the operation will be performed on.

        stat (string): the specific statistic to return. Possible statistics are "Count", "Unique Values", "First", "Last", "Sum", "Max", "Min", "Range", "Mean", "Median", "Mode", "IQ1", "IQ3", "IQR", "Standard Deviation", "Standard Error", "Variance", "Coefficient of Variation", "Skewness", and "Kurtosis"
        
        decimals (integer): the number of decimal places float values will be rounded to."""
        if not isinstance(row, int) or row < 0 or row >= self.num_rows():
            raise ValueError("Row index is out of range.")
        valid_stats = {
        "Count", "Unique Values", "First", "Last", "Sum", "Max", "Min", "Range", "Mean",
        "Median", "Mode", "IQ1", "IQ3", "IQR", "Standard Deviation",
        "Standard Error", "Variance", "Coefficient of Variation",
        "Skewness", "Kurtosis"}

        if stat not in valid_stats:
            raise ValueError(f"Unknown statistic '{stat}'.")

        r = self.get_row(row)
        r = _Math().remove_non_numeric(r)
        if not r:
            return None
        return _Stats().descriptive_stat(r,stat,decimals)
    
    def row_stat_list(self,stat,decimals = 2):
        """Calculate the specified descriptive statistic ('stat') of each row in the Table and return a list of statistics. Only numeric values (integer or float) in each row will be included in the calculation.

        Parameters:

        stat (string): the specific statistic to calculate for each row. Possible statistics are "Count", "Unique Values", "First", "Last", "Sum", "Max", "Min", "Range", "Mean", "Median", "Mode", "IQ1", "IQ3", "IQR", "Standard Deviation", "Standard Error", "Variance", "Coefficient of Variation", "Skewness", and "Kurtosis"
        
        decimals (integer): the number of decimal places float values will be rounded to."""
        valid_stats = {
        "Count", "Unique Values", "First", "Last", "Sum", "Max", "Min", "Range", "Mean",
        "Median", "Mode", "IQ1", "IQ3", "IQR", "Standard Deviation",
        "Standard Error", "Variance", "Coefficient of Variation",
        "Skewness", "Kurtosis"}

        if stat not in valid_stats:
            raise ValueError(f"Unknown statistic '{stat}'.")
        
        output = []
        stats = _Stats()
        maths = _Math()
        for row in self._data:
            r = maths.remove_non_numeric(row)
            if not r:
                output.append(None)
            else:
                output.append(stats.descriptive_stat(r,stat,decimals))
        return output
    
    def get_date_range(self,column):
        """Return a tuple of (start date, end date) from the specified column of dates.
        
        Parameters:
        
        column (integer, string): the index or header of the column the operation will be performed on."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in the table.")
        
        # Resolve column index
        if isinstance(column, str):
            column = self._headers.index(column)

        dates = self.get_column(column)

        return min(dates),max(dates)

    def date_to_iso(self,column,formats,new_header=None):
        """Convert a specified column ('column') to ISO format (yyyy-mm-dd) and either append a new column or overwrite an existing column.
        
        Parameters:
        
        column (integer, string): the index or header of the column the operation will be performed on.
        
        formats (string, list): a string or list of strings of date formats to attempt to convert to ISO format.
        
        new_header (string): if a string is given, a new column of converted ISO dates will be appended to the end of the data with the specified header. If left as None, the target date column will be overwritten with converted ISO dates."""
        # Validate column existence
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in the table.")

        # Resolve column index
        if isinstance(column, str):
            column = self._headers.index(column)

        # Normalize formats into a list
        if isinstance(formats, str):
            formats = [formats]
        elif not isinstance(formats, list):
            raise TypeError("formats must be a string or a list of strings.")
        
        strptime = datetime.datetime.strptime
        iso_fmt = "%Y-%m-%d"


        # Prepare list to store converted ISO dates
        iso_values = []

        for row in self.get_data():
            val = row[column]

            # Default: preserve original value
            iso_val = val

            if isinstance(val, str):
                for fmt in formats:
                    try:
                        dt = strptime(val, fmt)
                        iso_val = dt.strftime(iso_fmt)
                        break
                    except Exception:
                        continue

            iso_values.append(iso_val)


        # Overwrite mode
        if new_header is None:
            for i, row in enumerate(self._data):
                row[column] = iso_values[i]

        # Append mode — use built-in append_column
        else:
            self.append_column(iso_values, new_header)
        return self

    def normalize_dates(self, column, sep="-", new_header=None):
        """Normalize separators in a column of dates ('column') to a consistent ISO-style format and either append a new column or overwrite an existing column.

        Parameters:

        column (integer, string): the index or header of the column the operation will be performed on.

        sep (string): the separator to use in the normalized output. Default is '-'.
        
        new_header (string): if a string is given, a new column of normalized ISO dates will be appended to the end of the data with the specified header. If left as None, the target date column will be overwritten with normalized ISO dates."""
        dates = self.get_column(column)

        output = []
        for d in dates:
            # detect separator by finding the first non-digit character
            sep_in = next((ch for ch in d if not ch.isdigit()), None)

            if sep_in is None:
                # handle raw digit format 'yyyymmdd'
                if len(d) == 8 and d.isdigit():
                    year, month, day = d[:4], d[4:6], d[6:]
                    output.append(f"{year}{sep}{month}{sep}{day}")
                else:
                    raise ValueError(f"Invalid date format '{d}'. Expected 'yyyy{sep}mm{sep}dd'.")
            else:
                parts = d.split(sep_in)
                if len(parts) != 3:
                    raise ValueError(f"Invalid date format '{d}'. Expected three parts separated by '{sep_in}'.")
                year, month, day = parts
                output.append(f"{year}{sep}{month}{sep}{day}")

        # Overwrite mode
        if new_header is None:
            self.replace_column(column,output,new_header)

        # Append mode — use built-in append_column
        else:
            self.append_column(output, new_header)
        return self
    
    def check_date_continuity(self, column):
        """Check whether a column of dates ('column') forms a continuous, gap-free daily sequence from the earliest date to the latest date and return True or False.

        Parameters:

        column (integer, string): the index or header of the date column to be checked."""
        # Extract column as list
        dates = self.get_column(column)

        # Sort dates and get range
        sorted_dates = sorted(dates)
        start, end = self.get_date_range(sorted_dates)

        # generate_day_range is non-inclusive of end, so add one day manually
        end_dt = datetime.datetime.strptime(end, "%Y-%m-%d") + datetime.timedelta(days=1)
        end_plus = end_dt.strftime("%Y-%m-%d")

        full_range = _Date().generate_day_range(start, end_plus)

        return len(full_range) == len(sorted_dates)

    def insert_missing_date_values(self, column, placeholder=""):
        """Insert rows into the Table wherever dates are missing from the specified date column ('date_column'). The date column will receive the missing date, and all other columns will receive the specified placeholder.

        Parameters:

        column (integer, string): the index or header of the date column.

        placeholder: the value to insert for non-date columns in the new rows. Default is ""."""
        # Resolve column index
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")

        if isinstance(column, str):
            date_idx = self._headers.index(column)
        else:
            date_idx = column

        # Sort rows by date
        sorted_rows = sorted(self._data, key=lambda r: r[date_idx])
        sorted_dates = [row[date_idx] for row in sorted_rows]

        # Determine full expected date range
        start, end = _Date().get_date_range(sorted_dates)

        # generate_day_range is non-inclusive of end, so add one day manually
        end_dt = datetime.datetime.strptime(end, "%Y-%m-%d") + datetime.timedelta(days=1)
        end_plus = end_dt.strftime("%Y-%m-%d")

        full_range = _Date().generate_day_range(start, end_plus)

        # Build a map of existing rows by date
        row_map = {row[date_idx]: row for row in sorted_rows}
        date_set = set(sorted_dates)

        # Build new rows list
        new_rows = []
        lh = len(self._headers)

        for d in full_range:
            if d in date_set:
                new_rows.append(row_map[d])
            else:
                # missing date → create placeholder row
                row = [placeholder] * lh
                row[date_idx] = d
                new_rows.append(row)

        # build output table
        _out = TableTools().new_table()
        _out.set_data(new_rows)
        _out.set_headers(self.get_headers())
        _out.set_dtypes(self.get_dtypes())

        return _out

    # ============================================================
    # FILLING/INTERPOLATION
    # ============================================================
        
    def fill_data_single_value(self,column,fill_val,to_fill = ""):
        """Interpolate a new value for each missing or invalid value ('to_fill') in the specified column ('column') using a single specified value ('fill_val').

        Parameters:

        column (integer, string): the index or header of the column the operation will be performed on.

        fill_val (integer, float, string): the new value that will replace missing or invalid values.

        to_fill (integer, float, string, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in the table.")
      
        if isinstance(column,str):
            column = self._headers.index(column)
            
        if not isinstance(to_fill,list):
            to_fill = [to_fill]
        
        for row in self._data:
            if row[column] in to_fill:
                row[column] = fill_val
        
        self._dtypes[column] = self._return_col_type(self.get_column(column))
        return self
    
    def fill_data_from_column(self,column,other_col,to_fill = ""):
        """Interpolate a new value for each missing or invalid value ('to_fill') in the specified column ('column') using the corresponding row value from another column ('other_col'). The column data type will be redetermined.

        Parameters:

        column (integer, string): the index or header of the column the operation will be performed on.

        other_col (integer, string): the index or header of the column that will be used to fill missing values in the column of operation.

        to_fill (integer, float, string, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in the table.")
        if not self.column_exists(other_col):
            raise ValueError(f"Column '{other_col}' does not exist in the table.")
        
        if isinstance(column,str):
            column = self._headers.index(column)

        if isinstance(other_col,str):
            other_col = self._headers.index(other_col)
        
        if not isinstance(to_fill,list):
            to_fill = [to_fill]

        for row in self._data:
            if row[column] in to_fill and row[other_col] not in to_fill:
                row[column] = row[other_col]
        
        # we have to redetect the column data type in case there is a new value type
        self._dtypes[column] = self._return_col_type(self.get_column(column))
        return self

    def fill_data_from_list(self,column,fill_list,to_fill = ""):
        """Interpolate a new value for each missing or invalid value ('to_fill') in the specified column ('column') using the corresponding value from an input list ('fill_list'). The input list must be the same length as the column to be filled. The column data type will be redetermined.

        Parameters:

        column (integer, string): the index or header of the column the operation will be performed on.

        fill_list (list): the list of values that will be used to fill missing values in the column of operation. The input list must be the same length as the column to be filled.

        to_fill (integer, float, string, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in the table.")
        if not isinstance(fill_list,list):
            raise ValueError("Fill list must be a list.")
        if len(fill_list) != self.num_rows():
            raise ValueError("Length of input fill list must be equal to length of column.")
        
        if isinstance(column,str):
            column = self._headers.index(column)
        
        if not isinstance(to_fill,list):
            to_fill = [to_fill]

        for ix, row in enumerate(self._data):
            if row[column] in to_fill and fill_list[ix] not in to_fill:
                row[column] = fill_list[ix]
        
        # we have to redetect the column data type in case there is a new value type
        self._dtypes[column] = self._return_col_type(self.get_column(column))
        return self

    def fill_data_column_stat(self,column,stat,to_fill = "",decimals = 2):
        """Interpolate a new value for each missing or invalid value ('to_fill') in the specified column ('column') using the specified statistics ('stat') of all valid values within that column. Missing or invalid values will be interpolated as floats.

        Parameters:

        column (integer, string): the index or header of the column the operation will be performed on.

        stat (string): the specific statistic to return. Possible statistics are "Count", "Unique Values", "First", "Last", "Sum", "Max", "Min", "Range", "Mean", "Median", "Mode", "IQ1", "IQ3", "IQR", "Standard Deviation", "Standard Error", "Variance", "Coefficient of Variation", "Skewness", and "Kurtosis"
        
        to_fill (integer, float, string, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value.
        
        decimals (integer): the number of decimal places float values will be rounded to."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")
        
        valid_stats = {
        "Count", "Unique Values", "First", "Last", "Sum", "Max", "Min", "Range", "Mean",
        "Median", "Mode", "IQ1", "IQ3", "IQR", "Standard Deviation",
        "Standard Error", "Variance", "Coefficient of Variation",
        "Skewness", "Kurtosis"}

        if stat not in valid_stats:
            raise ValueError(f"Unknown statistic '{stat}'.")
        
        if isinstance(column,str):
            column = self._headers.index(column)
            
        if not isinstance(to_fill,list):
            to_fill = [to_fill]

        c = self.get_column(column)
        c = _Math().remove_non_numeric(c)
        if not c:
            raise ValueError(f"Column '{column}' has no valid values to compute statistic from.")
        stat = _Stats().descriptive_stat(c,stat,decimals)

        for row in self._data:
            if row[column] in to_fill:
                row[column] = stat
        
        self._dtypes[column] = self._return_col_type(self.get_column(column))
        return self
    
    def fill_data_forward(self,column,to_fill = ""):
        """Interpolate a new value for each missing or invalid value ('to_fill') in the specified column ('column') using the valid value immediately preceding each missing or invalid value. This method will not be able to fill missing or invalid values in the first row.

        Parameters:

        column (integer, string): the index or header of the column the operation will be performed on.
        
        to_fill (integer, float, string, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")
        if isinstance(column,str):
            column = self._headers.index(column)
        
        if not isinstance(to_fill,list):
            to_fill = [to_fill]
        
        for row in range(1, len(self._data)):
            if self._data[row][column] in to_fill:
                prev_val = self._data[row-1][column]
                if prev_val not in to_fill:
                    self._data[row][column] = prev_val
        return self
    
    def fill_data_backward(self,column,to_fill = ""):
        """Interpolate a new value for each missing or invalid value ('to_fill') in the specified column ('column') using the valid value immediately following each missing or invalid value. This method will not be able to fill missing or invalid values in the last row.

        Parameters:

        column (integer, string): the index or header of the column the operation will be performed on.
        
        to_fill (integer, float, string, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")
        if isinstance(column,str):
            column = self._headers.index(column)
        
        if not isinstance(to_fill,list):
            to_fill = [to_fill]
            
        for row in range(len(self._data) - 1):
            if self._data[row][column] in to_fill:
                next_val = self._data[row+1][column]
                if next_val not in to_fill:
                    self._data[row][column] = next_val 
        return self
    
    def fill_data_linear_interpolation(self,column,to_fill = "",direction = "forward"):
        """Interpolate a new value for each missing or invalid value ('to_fill') in the specified column ('column') using a linear interpolation technique. This method will not be able to fill missing or invalid values at either the beginning or end of the column, depending on the direction of interpolation.

        Parameters:

        column (integer, string): the index or header of the column the operation will be performed on.
        
        to_fill (integer, float, string, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value.
        
        direction (string): the direction of interpolation applied after the linear interpolation. May be specified as "forward" for a forward interpolation (missing or invalid values at the beginning of the column will not be interpolated) or "backward" for a backward interpolation (missing or invalid values at the end of the column will not be interpolated)."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")
        if isinstance(column,str):
            column = self._headers.index(column)
        
        if not isinstance(to_fill,list):
            to_fill = [to_fill]
        rows = self.num_rows()
        x_col = list(range(1, rows + 1))

        # Perform interpolation
        for i in range(rows):
            if self._data[i][column] in to_fill and 0 < i < rows - 1:
                # Find next valid value
                y2, x2 = None, None
                for j in range(i + 1, rows):
                    try:
                        y2 = float(self._data[j][column])
                        x2 = x_col[j]
                        break
                    except (ValueError, TypeError):
                        continue

                if y2 is not None:
                    try:
                        x, x1 = x_col[i], x_col[i - 1]
                        y1 = float(self._data[i - 1][column])
                        y = y1 + ((x - x1) / (x2 - x1)) * (y2 - y1)
                        self._data[i][column] = y
                    except (ValueError, TypeError, ZeroDivisionError):
                        pass

        if direction == "forward":
            self.fill_data_forward(column,to_fill)
        elif direction == "backward":
            self.fill_data_backward(column,to_fill)
        return self

    def fill_data_mean_distance(self,column,distance,to_fill = "",decimals = 2):
        """Interpolate a new value for each missing or invalid value ('to_fill') in the specified column ('column') using the average of valid values before and after each missing or invalid value. Missing or invalid values will be interpolated as floats. This technique is best suited to ordered data where the order is meaningful (e.g. data ordered by date or time).

        Parameters:

        column (integer, string): the index or header of the column the operation will be performed on.

        distance (integer): the number of values before and after each missing or invalid value to use for interpolation. Only valid (numeric) values within this range will be used for interpolation.

        to_fill (integer, float, string, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value.
        
        decimals (integer): the number of decimal places float values will be rounded to."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")
        if isinstance(column,str):
            column = self._headers.index(column)
        
        if not isinstance(to_fill,list):
            to_fill = [to_fill]
        
        if distance <= 0:
            raise ValueError(f"Distance value must be a positive integer.")

        rows = len(self._data)

        for i in range(rows):
            if self._data[i][column] in to_fill:
                vals_for_mean = []
                for d in range(1, distance + 1):
                    # Look forward
                    if i + d < rows:
                        val = self._data[i + d][column]
                        if val not in to_fill:
                            try:
                                vals_for_mean.append(float(val))
                            except (ValueError, TypeError):
                                pass
                    # Look backward
                    if i - d >= 0:
                        val = self._data[i - d][column]
                        if val not in to_fill:
                            try:
                                vals_for_mean.append(float(val))
                            except (ValueError, TypeError):
                                pass

                if vals_for_mean:
                    self._data[i][column] = round(sum(vals_for_mean) / len(vals_for_mean), decimals)
        return self

    def fill_data_rolling_average(self, column, window, to_fill="", decimals=2):
        """Interpolate a new value for each missing or invalid value ('to_fill') in the specified column ('column') using the average of valid values within a rolling window of previous rows. Missing or invalid values will be interpolated as floats. This technique is best suited to ordered data where the order is meaningful (e.g. data ordered by date or time).

        Parameters:

        column (integer, string): the index or header of the column the operation will be performed on.

        window (integer): the number of previous values to use for interpolation. Only valid (numeric) values within this range will be used for interpolation.

        to_fill (integer, float, string, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value.
        
        decimals (integer): the number of decimal places float values will be rounded to."""
        # Validate column
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in the table.")

        # Resolve column index
        if isinstance(column, str):
            column = self._headers.index(column)

        # Normalize to_fill
        if not isinstance(to_fill, list):
            to_fill = [to_fill]

        # Validate window
        if window <= 0:
            raise ValueError("Window size must be a positive integer.")

        rows = len(self._data)

        # Perform rolling average fill
        for i in range(rows):
            if self._data[i][column] in to_fill:
                vals_for_mean = []
                for j in range(max(0, i - window), i):  # look only backward
                    val = self._data[j][column]
                    if val not in to_fill:
                        try:
                            vals_for_mean.append(float(val))
                        except (ValueError, TypeError):
                            pass
                if vals_for_mean:
                    self._data[i][column] = round(sum(vals_for_mean) / len(vals_for_mean), decimals)
        return self

    def fill_data_inverse_distance_weighted(self,column,distance,power = 2,to_fill = "",decimals = 2):
        """Interpolate a new value for each missing or invalid value ('to_fill') in the specified column ('column') using an inverse distance weighted (IDW) technique. Missing or invalid values will be interpolated as floats. This technique is best suited to ordered data where the order is meaningful (e.g. data ordered by date or time).

        Parameters:

        column (integer, string): the index or header of the column the operation will be performed on.

        distance (integer): the number of values before and after each missing or invalid value to use for interpolation. Only valid (numeric) values within this range will be used for interpolation.

        power (integer): the power of the inverse distance weighting function. Higher values will give greater weight to nearer values i.e. nearer values will contribute more strongly to the interpolated value.

        to_fill (integer, float, string, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value.
        
        decimals (integer): the number of decimal places float values will be rounded to."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in Table.")
        if isinstance(column,str):
            column = self._headers.index(column)
        
        if not isinstance(to_fill,list):
            to_fill = [to_fill]
        
        if distance <= 0:
            raise ValueError(f"Distance value must be a positive integer.")
        rows = len(self._data)

        # IDW equation
        # v1/(d1**power) + vn/(dn**power) / 1/(d1**power) + 1/(dn**power)

        for i in range(rows):
            if self._data[i][column] in to_fill:
                numerators = []
                denominators = []
                for d in range(1, distance + 1):
                    # Look forward
                    if i + d < rows:
                        val = self._data[i + d][column]
                        if val not in to_fill:
                            try:
                                numerators.append(float(val) / (d ** power))
                                denominators.append(1 / (d ** power))
                            except (ValueError, TypeError):
                                pass
                    # Look backward
                    if i - d >= 0:
                        val = self._data[i - d][column]
                        if val not in to_fill:
                            try:
                                numerators.append(float(val) / (d ** power))
                                denominators.append(1 / (d ** power))
                            except (ValueError, TypeError):
                                pass

                if numerators and denominators:
                    self._data[i][column] = round(sum(numerators) / sum(denominators), decimals)
        return self
    
    def fill_data_regression(self,column,vals,to_fill = "",decimals = 2):
        """Interpolate a new value for each missing or invalid value ('to_fill') in the specified column ('column') by performing a linear regression between the values of the specified column and a second list ('vals'). Missing or invalid values will be interpolated as floats.

        Parameters:

        column (integer, string): the index or header of the column the operation will be performed on.

        vals (list): the list of values to perform a linear regression with for interpolation. Only valid (numeric) values within this list will be used for interpolation. Should be the same length or longer than the specified column.

        to_fill (integer, float, string, list): the value or list of values to be considered missing or invalid. The default is an empty string "" representing a missing value.
        
        decimals (integer): the number of decimal places float values will be rounded to."""
        if not self.column_exists(column):
            raise ValueError(f"Column '{column}' does not exist in the table.")
        if not isinstance(vals,list):
            raise ValueError("Fill list must be a list.")
        if len(vals) != self.num_rows():
            raise ValueError("Length of input fill list must be equal to length of column.")
        
        y_vals, x_vals = [], []
        for i in range(len(self._data)):
            val = self._data[i][column]
            if val not in to_fill:
                try:
                    y_vals.append(float(val))
                    x_vals.append(float(vals[i]))
                except (ValueError, TypeError):
                    pass

        if not y_vals or not x_vals:
            raise ValueError("No valid numeric values available for regression.")

        # Calculate regression coefficients
        model = _Stats().linear_regression(x_vals, y_vals)
        intercept, slope = model["intercept"],model["slope"]

        # Fill gaps using regression prediction
        for i in range(len(self._data)):
            if self._data[i][column] in to_fill and vals[i] not in to_fill:
                try:
                    self._data[i][column] = round((float(vals[i]) * slope) + intercept, decimals)
                except (ValueError, TypeError):
                    pass
        return self

    # ============================================================
    # DIRECT PLOTTING
    # ============================================================

    def scatter_plot(self,x_col, y_cols,title="Scatter Plot",x_label="Independent Variable",y_label="Dependent Variable",num_ticks=5,point_size=3,point_colour=None, point_style=None,trendline=False,trend_colour=None,legend=False,point_labels=None,tooltip_threshold=100,margin=60, filename="scatter_plot.html",open_browser=False):
        """Create a scatter plot for a single independent variable ('x_col') column and one or more dependent variables ('y_cols') and configure plot elements.

        Parameters:

        x_col (integer, string): the index or header of the column representing the independent variable to be visualized. Values may be numerical (integer or floating point) or categorical strings such as dates.

        y_cols (integer, string, list): the index or header of a single column representing the dependent variable to be visualized, or a list of columns indices or headers if there are multiple dependent variables to be visualized. Values must be numerical (integer or floating point).

        title (string): the title of the scatter plot. 

        x_label (string): the label of the x axis.

        y_label (string): the label of the y axis. 

        num_ticks (integer): the number of tick marks to draw on each axis.

        point_size (integer): the radius of points on the plot (in pixels).

        point_colour (string, list): if a single dependent variable is given, this should be a single string representing the colour of the points. If multiple dependent variables are given, this should be a list of strings corresponding to the colour of each variable. If no colours are specified or the number of specified colours does not equal the number of dependent variables, random colours will be generated. Colours may be specified as a common colour name or hex code.

        point_style (string, list): the shape of points. May be specified as "circle", "diamond", "square", or "triangle". If a single dependent variable is given, this should be a single string representing the style of the points. If multiple dependent variables are given, this should be a list of strings corresponding to the style of each variable. If no styles are specified or the number of specified styles does not equal the number of dependent variables, random styles will be generated.

        trendline (Boolean): flag to indicate if trend lines for each dependent variable will be calculated and displayed.

        trend_colour (string, list): if a single dependent variable is given, this should be a single string representing the colour of the trend line. If multiple dependent variables are given, this should be a list of strings corresponding to the colour of each trend line. If no colours are specified or the number of specified colours does not equal the number of dependent variables, random colours will be generated. Colours may be specified as a common colour name or hex code.
        
        legend (Boolean): flag to indicate if a legend will be displayed.

        point_labels (string, list): if a single dependent variable is given, this should be a single string representing the label of the points. If multiple dependent variables are given, this should be a list of strings representing the labels of each variable. If no labels are given, generic labels will be generated. A legend must be created to view point labels.

        tooltip_threshold (integer): the maximum number of points per series for which explicit hover-based tooltips are included.
        
        margin (integer): the margin size (in pixels) around the plot area. 
        
        filename (string): the name of the output file. Include the output file path in the string or prepend tt.get_out_dir() to the filename.

        open_browser (Boolean): flag to indicate if the generated file should be automatically opened in the default web browser."""
        x_col = self.get_column(x_col)
        if isinstance(y_cols,(int,str)):
            y_cols = self.get_column(y_cols)
            if not all(isinstance(v, (int, float)) for v in y_cols):
                raise ValueError("All values in 'y_cols' must be numeric (int or float).")
        elif isinstance(y_cols,list):
            y_cols = self.get_columns(y_cols)
            for series in y_cols:
                if not all(isinstance(v, (int, float)) for v in series):
                    raise ValueError("All values in each dependent variable series must be numeric (int or float).")  
        _Plot().scatter_plot(x_col,y_cols,title,x_label,y_label,num_ticks,point_size,point_colour,point_style,trendline,trend_colour,legend,point_labels,tooltip_threshold,margin,filename,open_browser)
        return self

    def line_plot(self,x_col, y_cols,title="Line Plot",x_label="Independent Variable", y_label="Dependent Variable",num_ticks=5,line_weight=2,line_colours=None,trendline=False,trend_colours=None,legend=False,line_labels=None,margin=60,filename="line_plot.html",open_browser=False):
        """Create a line plot for a single independent variable ('x_col') column and one or more dependent variables ('y_col') and configure plot elements.

        Parameters:

        x_col (integer, string): the index or header of the column representing the independent variable to be visualized. Values may be numerical (integer or floating point) or categorical strings such as dates.

        y_cols (integer, string, list): the index or header of a single column representing the dependent variable to be visualized, or a list of columns indices or headers if there are multiple dependent variables to be visualized. Values must be numerical (integer or floating point).

        title (string): the title of the line plot.  

        x_label (string): the label of the x axis.

        y_label (string): the label of the y axis.  

        num_ticks (integer): the number of tick marks to draw on each axis.

        line_weight (integer): the stroke width of the lines on the plot (in pixels).

        line_colours (string, list): if a single dependent variable is given, this should be a single string representing the colour of the line. If multiple dependent variables are given, this should be a list of strings corresponding to the colour of each variable. If no colours are specified or the number of specified colours does not equal the number of dependent variables, random colours will be generated. Colours may be specified as a common colour name or hex code.

        trendline (Boolean): flag to indicate if trend lines for each dependent variable will be calculated and displayed.

        trend_colours (string, list): if a single dependent variable is given, this should be a single string representing the colour of the trend line. If multiple dependent variables are given, this should be a list of strings corresponding to the colour of each trend line. If no colours are specified or the number of specified colours does not equal the number of dependent variables, random colours will be generated. Colours may be specified as a common colour name or hex code.
        
        legend (Boolean): flag to indicate if a legend will be displayed.

        line_labels (string, list): the labels of each dependent variable series. If no labels are given, generic labels ("Series 1", "Series 2", ...) will be generated. A legend must be created to view series labels.
        
        margin (integer): the margin size (in pixels) around the plot area.  
        
        filename (string): the name of the output file. Include the output file path in the string or prepend tt.get_out_dir() to the filename.

        open_browser (Boolean): flag to indicate if the generated file should be automatically opened in the default web browser."""
        x_col = self.get_column(x_col)
        if isinstance(y_cols,(int,str)):
            y_cols = self.get_column(y_cols)
            if not all(isinstance(v, (int, float)) for v in y_cols):
                raise ValueError("All values in 'y_cols' must be numeric (int or float).")
        elif isinstance(y_cols,list):
            y_cols = self.get_columns(y_cols)
            for series in y_cols:
                if not all(isinstance(v, (int, float)) for v in series):
                    raise ValueError("All values in each dependent variable series must be numeric (int or float).") 
        _Plot().line_plot(x_col,y_cols,title,x_label,y_label,num_ticks,line_weight,line_colours,trendline,trend_colours,legend,line_labels,margin,filename,open_browser)
        return self

    def histogram(self,column,title="Histogram",x_label="Value",y_label="Frequency",num_bins=None,bar_colour=None,margin=60,filename="histogram.html",open_browser=False):
        """Create a histogram for a single column and configure plot elements.

        Parameters:

        column (integer, string): the index or header of the column to be visualized.

        title (string): the title of the histogram.

        x_label (string): the label of the x axis.

        y_label (string): the label of the y axis.

        num_bins (integer, optional): the number of bins to divide the data into. If not specified, the optimal number of bins is calculated using Sturge's Formula.

        bar_colour (string, optional): the colour of the bars. If not specified, a random colour will be generated. Colours may be specified as a common colour name or hex code.

        margin (integer): the margin size (in pixels) around the plot area.

        filename (string): the name of the output file. Include the output file path in the string or prepend tt.get_out_dir() to the filename.

        open_browser (Boolean): flag to indicate if the generated file should be automatically opened in the default web browser."""
        column = self.get_column(column)
        if not all(isinstance(v, (int, float)) for v in column):
            raise ValueError("All values in 'column' must be numeric (int or float).")
        
        _Plot().histogram(column,title,x_label,y_label,num_bins,bar_colour,margin,filename,open_browser)
        return self
       
    def box_plot(self,category_cols,labels=None,title="Box Plot",x_label="Categories",y_label="Values",num_ticks=5,colour="#4e79a7",margin=60,filename="box_plot.html",open_browser=False):
        """Create a box plot from a list of columns ('category_cols') and configure plot elements. One box will be created per column.

        Parameters:

        category_cols (list): the nested list of column indices or headers of each category to be visualized.

        labels (list): the list of category labels for the x axis. Must be the same length as category_cols.

        title (string): the title of the box plot.

        x_label (string): the label of the x axis.

        y_label (string): the label of the y axis.

        num_ticks (integer): the number of tick marks to draw on the y axis.

        colour (string): the colour of the boxes. If not specified, a random colour will be generated. Colours may be specified as a common colour name or hex code.

        margin (integer): the margin size (in pixels) around the plot area.

        filename (string): the name of the output file. Include the output file path in the string or prepend tt.get_out_dir() to the filename.

        open_browser (Boolean): flag to indicate if the generated file should be automatically opened in the default web browser."""
        if isinstance(category_cols,(int,str)):
            category_cols = [category_cols]
        if isinstance(labels,str):
            labels = [labels]
        columns = self.get_columns(category_cols)
        
        _Plot().box_plot(columns,labels,title,x_label,y_label,num_ticks,colour,margin,filename,open_browser)
        return self
    
    def bar_chart(self,column,title="Bar Chart",x_label="Category",y_label="Count",num_ticks=6,colour="#4e79a7",margin=60,filename="bar_chart.html",open_browser=False):
        """Create a bar chart from a column of categorical values ('column') and configure chart elements.

        Parameters:

        column (integer, string): the index or header of the column to be visualized. 

        title (string): the title of the bar chart.

        x_label (string): the label of the x axis.

        y_label (string): the label of the y axis.

        num_ticks (integer): The number of tick marks to draw on the y axis.

        colour (string): the colour of the bars. If not specified, a random colour will be generated. Colours may be specified as a common colour name or hex code.

        margin (integer): the margin size (in pixels) around the plot area.

        filename (string): the name of the output file. Include the output file path in the string or prepend tt.get_out_dir() to the filename.

        open_browser (Boolean): flag to indicate if the generated file should be automatically opened in the default web browser."""
        column = self.get_column(column)
        _Plot().bar_chart(column,title,x_label,y_label,num_ticks,colour,margin,filename,open_browser)
        return self
    
    def pie_chart(self,column,title="Pie Chart",show_values=False,colours=None,legend=True,margin=60,filename="pie_chart.html",open_browser=False):
        """Create a pie chart from a column of categorical values ('column') and configure chart elements.

        Parameters:

        column (integer, string): the index or header of the column to be visualized. 

        title (string): the title of the pie chart.

        show_values (Boolean): flag to indicate if counts and percentages should be shown in the legend beside each category label. If False, only category names are shown.

        colours (list): list of colour strings to use for slices. If None, colours are auto-generated. Colours may be specified as a common colour name or hex code.

        legend (Boolean): flag to indicate if a legend will be displayed.

        margin (integer): the margin size (in pixels) around the plot area.

        filename (string): the name of the output file. Include the output file path in the string or prepend tt.get_out_dir() to the filename.

        open_browser (Boolean): flag to indicate if the generated file should be automatically opened in the default web browser."""
        column = self.get_column(column)
        _Plot().pie_chart(column,title,show_values,colours,legend,margin,filename,open_browser)
        return self
    
# TableTools class
class TableTools():
    """An object for accessing the functionality of the TableTools object and its toolboxes."""
    # ============================================================
    # PRIVATE METHODS AND CLASS ATTRIBUTES
    # ============================================================
    def __init__(self):
        self._help_index = self._load_help_index()
        self._in_dir = ""
        self._out_dir = ""
        self._start_time = None
        self._prog_count = 0
        self._prog_total = None
        self._curr_percent = 0
        self._prog_start_time = None
        self._allow_overwrite = False
        self._boolean_true_values = {v.lower() for v in ["true", "t", "yes"]}
        self._boolean_false_values = {v.lower() for v in ["false", "f", "no"]}
        self._all_extensions = [
        ".txt", ".dat", ".prn", ".fw", ".asc",
        ".csv", ".tsv", ".psv",
        ".json",
        ".html", ".htm",
        ".xml",
        ".yaml", ".yml",
        ".ini", ".cfg", ".conf",
        ".md", ".markdown",
        ".tex",
        ".r", ".dput",
        ".dbf",
        ".db", ".db2", ".db3",
        ".sqlite", ".sqlite2", ".sqlite3",]
        self._total_funcs = len([f for f in dir(__class__) if not f.startswith('_')])
        self.listops = _ListOps()
        self.generate = _Generate()
        self.math = _Math()
        self.stats = _Stats()
        self.conversion = _Conversion()
        self.points = _Points()
        self.text = _Text()
        self.date = _Date()
        self.plot = _Plot()
        self.matrix = _Matrix()
        self.symbols = _Symbols()
        self._all_toolboxes = [self.listops,self.generate,self.math,self.stats,self.conversion,self.points,self.text,self.date,self.plot,self.matrix]

    def __repr__(self):
        """If the TableTools object is returned, display a message."""
        return f"Object for accessing functionality of TableTools and its toolboxes."

    def _load_help_index(self):
        """Attempt to import the TableTools help index file and store it."""
        import json
        package_dir = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(package_dir, "tabletools_help_index.json")
        try:
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            return {}

    def _lookup_help_entry(self, topic):
        """Return the help entry dict for a given topic, or None if not found.
        
        Parameters:
        
        topic (string): the topic to search for and return."""
        if not topic:
            return None

        topic = topic.strip()

        # Exact match
        if topic in self._help_index:
            return self._help_index[topic]

        # Case-insensitive match
        lowered = {k.lower(): k for k in self._help_index}
        if topic.lower() in lowered:
            return self._help_index[lowered[topic.lower()]]

        # Prefix match (e.g., "math" is a toolbox)
        matches = [k for k in self._help_index if k.startswith(topic)]
        if len(matches) == 1:
            return self._help_index[matches[0]]

        # Fuzzy match (contains)
        contains = [k for k in self._help_index if topic.lower() in k.lower()]
        if len(contains) == 1:
            return self._help_index[contains[0]]

        # No match
        return None

    def _print_help_entry(self, entry):
        """Pretty-print a help entry.
        
        Parameters:
        
        entry (dictionary): the dictionary retrieved from the help index."""
        if not entry:
            print("No help available.")
            return

        name = entry.get("name", "")
        category = entry.get("category", "")
        signature = entry.get("signature", "")
        doc = entry.get("doc", "")
        example = entry.get("example", "")
        functions = entry.get("functions", None)

        print(f"\n=== {name} ===")

        # Human-friendly source labeling
        if category == "TableTools":
            print("Source: TableTools")
        elif category == "Table":
            print("Source: Table")
        elif category not in ("", "category"):
            print(f"Toolbox: {category}")

        if signature:
            print(f"Signature: {signature}")

        if doc:
            print("\nDescription:")
            print(doc)

        # Toolbox function list
        if functions:
            print("\nFunctions:")
            for fn in functions:
                print(f"    {fn}")

        if example:
            print("\nExample:")
            print(f"    {example}")

        print()

    def _row_dtype(self, row):
        """Convert the data type of each value in row when reading in data.
        
        Parameters:
        
        row (list): the row to convert raw string values into Python types."""
        
        out_row = []

        for r in row:
            # Strip whitespace and surrounding quotes
            if isinstance(r, str):
                r = r.strip()
                if r.startswith('"') and r.endswith('"'):
                    r = r[1:-1]

            # Boolean detection (case-insensitive)
            if isinstance(r, str):
                lower = r.lower()
                if lower in self._boolean_true_values:
                    out_row.append(True)
                    continue
                if lower in self._boolean_false_values:
                    out_row.append(False)
                    continue

            # Integer detection
            try:
                out_row.append(int(r))
                continue
            except (ValueError, TypeError):
                pass

            # Float detection
            try:
                out_row.append(float(r))
                continue
            except (ValueError, TypeError):
                pass

            # Fallback: string
            out_row.append(str(r))

        return out_row

    def _filter_in_dir(self, file_ext="", include="", exclude=""):
        """Return a list of files in the input directory after applying filters.
        
        Parameters:

        file_ext (string, list): the extension of the file type to be returned, or a list of extensions to be returned. May be left blank to include all file extensions.
        
        include (string, list): a character, subtring, or list of characters and/or substrings that must be included in the file name for the file to be returned. May be left blank to leave files unfiltered.
        
        exclude (string, list): a character, subtring, or list of characters and/or substrings that must not be in the file name for the file to be returned. May be left blank to leave files unfiltered."""
        # normalize to lists
        if isinstance(file_ext, str): file_ext = [file_ext] if file_ext else []
        if isinstance(include, str): include = [include] if include else []
        if isinstance(exclude, str): exclude = [exclude] if exclude else []

        files = []
        for f in os.listdir(self._in_dir):
            # extension filter
            if file_ext and not any(f.endswith(ext) for ext in file_ext):
                continue
            # include filter
            if include and not any(i in f for i in include):
                continue
            # exclude filter
            if exclude and any(e in f for e in exclude):
                continue
            files.append(f)
        return files

    def _format_time(self,seconds, decimals):
            """Format time in seconds and return the formatted time.
            
            Parameters:
            
            seconds (integer): the time to be formatted.
            
            decimals (integer): the number of decimal places time values will be rounded to."""
            units = "seconds"
            if seconds / 3600 > 1:
                seconds /= 3600
                units = "hours"
            elif seconds / 60 > 1:
                seconds /= 60
                units = "minutes"
            if seconds == 1:
                units = units[:-1]
            return f"{round(seconds, decimals)} {units}"

    def _extension_handler(self, filename, target_ext, valid_exts=None):
        """Ensure 'filename' ends with a valid extension.
        
        Parameters:
        
        filename (string): the name of the file to ensure the extension of.
        
        target_ext (string): the extension to enforce.
        
        valid_exts (list, string): other valid extensions for this file type."""
        lower = filename.lower()
        # Normalize valid_exts
        if valid_exts is None:
            valid_exts = [target_ext]
        elif isinstance(valid_exts, str):
            valid_exts = [valid_exts]

        # 1. Already valid then return unchanged
        for ext in valid_exts:
            if lower.endswith(ext.lower()):
                return filename

        # 2. Replace any known extension
        for ext in self._all_extensions:
            if lower.endswith(ext.lower()):
                return filename[: -len(ext)] + target_ext

        # 3. No extension then append target_ext
        return filename + target_ext

    def _read_logical_rows(self, file_obj):
        """Yield complete logical rows, assembling multi-line fields when needed.
        
        Parameters:
        
        file_object (object): the object returned by the file opener."""
        buffer = ""

        for line in file_obj:
            # Strip only the newline; preserve everything else
            line = line.rstrip("\n")

            if not buffer:
                buffer = line
            else:
                buffer += "\n" + line

            # Count quotes to determine if row is complete
            # Count *unescaped* quotes only
            quote_count = 0
            i = 0
            length = len(buffer)
            while i < length:
                if buffer[i] == '"':
                    # Escaped quote ("") → skip the next quote
                    if i + 1 < length and buffer[i + 1] == '"':
                        i += 1
                    else:
                        quote_count += 1
                i += 1

            # If quote_count is even → row is complete
            if quote_count % 2 == 0:
                yield buffer
                buffer = ""

        # If file ends while still in an incomplete row, yield what we have
        if buffer:
            yield buffer

    def _split_quoted_line(self, line, delimiter=","):
        """Split a delimited line into fields, respecting quoted fields.

        Parameters:

        line (list): the line of the input file to split.

        delimiter (string): the delimiter to split on."""
        fields = []
        current = []
        in_quotes = False
        i = 0
        length = len(line)

        while i < length:
            ch = line[i]

            if ch == '"':
                # Toggle quote state unless this is an escaped quote ("")
                if in_quotes and i + 1 < length and line[i + 1] == '"':
                    # Escaped quote: add a literal quote and skip the next char
                    current.append('"')
                    i += 1
                else:
                    # Enter or exit quoted mode
                    in_quotes = not in_quotes

            elif ch == delimiter and not in_quotes:
                # Delimiter outside quotes → end of field
                fields.append("".join(current))
                current = []

            else:
                # Normal character
                current.append(ch)

            i += 1

        # Append the final field
        fields.append("".join(current))

        return fields

    def _escape_field(self, value, delimiter):
        """Ensure a value is written correctly to a file if the value has commas.
        
        Parameters:
        
        value (integer, float, string, Boolean): the value to check.
        
        delimiter (string): the delimiter to escape."""
        s = str(value)

        # If it contains a quote, escape it by doubling
        if '"' in s:
            s = s.replace('"', '""')

        # If it contains comma, quote, or newline → wrap in quotes
        if delimiter in s or '"' in s or '\n' in s or '\r' in s:
            s = f'"{s}"'

        return s

    # ============================================================
    # CORE INFRASTRUCTURE & METADATA
    # ============================================================

    def version(self):
        """Print the current version of TableTools."""
        print(f"TableTools v1.0.4")

    def view_manual(self):
        """Open a web browser to view the TableTools manual."""
        import webbrowser as wb
        package_dir = os.path.dirname(os.path.abspath(__file__))
        docs_dir = os.path.join(package_dir, "Docs")
        manual_path = os.path.join(docs_dir, "TableToolsManual.html")
        wb.open_new_tab(manual_path)
    
    def help(self, topic=None):
        """Display the help entry for a given topic. The topic may refer to any function belonging to the TableTools object, a Table object, or any toolbox. If no topic is provided, a list of all available help topics will be printed.

        Parameters:

        topic (string): the name of the function or toolbox to display help for. Do not include function parentheses."""
        import difflib
        if topic is None:
            print("TableTools Help\n")
            print("Available topics:\n")
            for key in sorted(self._help_index):
                print(f"  {key}")
            print("\nUse tt.help('topic') for details.")
            return

        entry = self._lookup_help_entry(topic)

        if entry:
            self._print_help_entry(entry)
            return

        # No exact match then show suggestions
        print(f'No help available for "{topic}".\n')

        topic_lower = topic.lower()
        keys = list(self._help_index.keys())

        # 1. Substring matches
        substring_matches = [
            k for k in keys if topic_lower in k.lower()
        ]

        # 2. Prefix matches
        prefix_matches = [
            k for k in keys if k.lower().startswith(topic_lower)
        ]

        # 3. Close matches (difflib)
        close_matches = difflib.get_close_matches(topic, keys, n=5, cutoff=0.5)

        # Combine and deduplicate
        suggestions = []
        for group in (substring_matches, prefix_matches, close_matches):
            for k in group:
                if k not in suggestions:
                    suggestions.append(k)

        if suggestions:
            print("Did you mean:\n")
            for s in suggestions:
                print(f"    {s}")
            print()
        else:
            print("No similar topics found.\n")

    def display_total_functions(self):
        """Print the total number of functions in the TableTools object, the Table class, and all toolboxes."""
        total = self._total_funcs + _Table()._total_funcs
        for s in self._all_toolboxes:
            total += s._total_funcs
        print(f"Total number of functions in TableTools object, Table object, and toolboxes: {total}")
    
    def display_object_functions(self):
        """Print a list of all TableTools object methods."""
        print(f"Number of {__class__.__name__} functions: {len([f for f in dir(__class__) if not f.startswith('__') and not f.startswith('_')])}")
        for f in [f for f in dir(__class__) if not f.startswith("__") and not f.startswith('_')]:
            print(f)
    
    def display_toolboxes(self):
        """Print a list of all toolboxes accessed by the TableTools object, along with a general description of toolbox functionality."""
        for s in self._all_toolboxes:
            print(s)
   
    # ============================================================
    # CONFIGURATION AND STATE
    # ============================================================

    def set_overwrite_state(self,state):
        """Set the state for allowing TableTools to overwrite files when writing.
        
        Parameters:
        
        state (Boolean): flag to indicate if file overwrite is enabled."""
        self._allow_overwrite = state

    def get_overwrite_state(self):
        """Return the state for allowing TableTools to overwrite files when writing."""
        if self._allow_overwrite:
            return True
        else:
            return False
    
    def set_boolean_values(self, true_values=["true", "t", "yes"], false_values=["false", "f", "no"]):
        """Set the string values that TableTools will interpret as Boolean True or False. Values are case-insensitive.

        Parameters:

        true_values (list): list of strings to interpret as Boolean True.
        
        false_values (list): list of strings to interpret as Boolean False."""
        
        # Reset to defaults (lowercased)
        self._boolean_true_values = {v.lower() for v in ["true", "t", "yes"]}
        self._boolean_false_values = {v.lower() for v in ["false", "f", "no"]}

        # Add user-provided True values
        for v in true_values:
            if isinstance(v, str):
                self._boolean_true_values.add(v.lower())

        # Add user-provided False values
        for v in false_values:
            if isinstance(v, str):
                self._boolean_false_values.add(v.lower())

    def get_boolean_values(self):
        """Get the current string values that TableTools interprets as Boolean True or False. Returns a list of [boolean_true_values, boolean_false_values].""" 
        return [sorted(list(self._boolean_true_values)),sorted(list(self._boolean_false_values))]

    def set_in_dir(self,path):
        """Set an input directory ('path') to read files from. Setting an input directory is not required when reading files, however, if an input directory has not been set a path must be specified along with the filename when reading files.

        Parameters:
        
        path (string): the path of the input directory."""
        if not isinstance(path, str):
            raise TypeError("Input directory must be a string.")

        # Normalize empty or None, clear directory
        if not path:
            self._in_dir = ""
            return

        # Strip trailing separators and re-add exactly one
        path = path.rstrip(os.sep) + os.sep

        self._in_dir = path
    
    def get_in_dir(self):
        """Return the input directory if an input directory has been set."""
        d = self._in_dir
        return d

    def set_out_dir(self,path):
        """Set an output directory ('path') to write files to. Setting an output directory is not required when writing files, however, if an output directory has not been set a path must be specified along with the filename when writing files.

        Parameters:
        
        path (string): the path of the output directory."""
        if not isinstance(path, str):
            raise TypeError("Output directory must be a string.")

        # Normalize empty or None, clear directory
        if not path:
            self._out_dir = ""
            return

        # Strip trailing separators and re-add exactly one
        path = path.rstrip(os.sep) + os.sep

        self._out_dir = path
    
    def get_out_dir(self):
        """Return the output directory if an output directory has been set."""
        d = self._out_dir
        return d

    # ============================================================
    # DIRECTORY AND FILE MANAGEMENT
    # ============================================================

    def create_directory(self,dir_name):
        """Attempt to create a directory ('dir_name').
        
        Parameters:
        
        dir_name (string): the name of the directory to create. No directory will be created if a folder with this name already exists in this location."""
        if not isinstance(dir_name, str):
            raise TypeError("Directory path must be a string.")

        # Normalize path (remove trailing separators)
        dir_name = dir_name.rstrip(os.sep)

        try:
            os.makedirs(dir_name, exist_ok=False)
        except FileExistsError:
            print("A folder with this name already exists at this location.")

    def open_file(self,filename=""):
        """Open a file ('filename') with the default program specified by the operating system.
        
        Parameters:
        
        filename (string): the name of the file to be opened, including the path and extension."""
        import subprocess
        import sys

        if not isinstance(filename, str):
            raise TypeError("File path must be a string.")

        # Normalize to an absolute path
        filename = os.path.abspath(filename)

        # Cross‑platform file opening
        try:
            if sys.platform.startswith("win"):
                os.startfile(filename)
            elif sys.platform.startswith("darwin"):
                subprocess.call(["open", filename])
            else:
                subprocess.call(["xdg-open", filename])
        except Exception as e:
            print("Error opening file:", e)

    def file_exists(self,filename):
        """Determine if a file ('filename') exists. Returns True or False.

        Parameters:

        filename (string): the name of the file to check for, including the path and extension."""
        if not isinstance(filename, str):
            raise TypeError("File path must be a string.")

        # Normalize to an absolute path
        filename = os.path.abspath(filename)

        return os.path.exists(filename)
 
    def directory_exists(self,path):
        """Determine if a directory ('path') exists. Returns True or False.

        Parameters:

        path (string): the name of the directory to check for."""
        if not isinstance(path, str):
            raise TypeError("Directory path must be a string.")

        # Normalize to an absolute path
        path = os.path.abspath(path)

        return os.path.isdir(path)

    def delete_file(self,filename): 
        """Attempt to delete a file ('filename'). Files are NOT sent to the recycle bin and cannot be recovered once deleted.

        Parameters:

        filename (string): the name of the file to delete, including the path and extension."""
        if not isinstance(filename, str):
            raise TypeError("File path must be a string.")

        # Normalize to an absolute path
        filename = os.path.abspath(filename)

        if os.path.isfile(filename):
            try:
                os.remove(filename)
            except Exception as e:
                print(f"Error deleting file: {e}")
        else:
            print(f"File '{filename}' does not exist.")

    def delete_directory(self,path): 
        """Attempt to delete a directory ('path') and all its contents. Directory and contents are NOT sent to the recycle bin and cannot be recovered once deleted.

        Parameters:

        path (string): the name of the directory to delete for, including the path and extension."""
        import shutil

        if not isinstance(path, str):
            raise TypeError("Directory path must be a string.")

        # Normalize to an absolute path
        path = os.path.abspath(path)

        if os.path.isdir(path):
            try:
                shutil.rmtree(path)
            except Exception as e:
                print(f"Error deleting directory: {e}")
        else:
            print(f"Directory '{path}' does not exist.")
  
    def num_files_in_input_dir(self,file_ext = "",include = "",exclude = ""):
        """Return the number of files in the input directory.

        Parameters:

        file_ext (string, list): the extension of the file type to be counted. May be left blank to include all file extensions.
        
        include (string, list): a character, subtring, or list of characters and/or substrings that must be included in the file name for the file to be counted. May be left blank to leave files unfiltered..
        
        exclude (string, list): a character, subtring, or list of characters and/or substrings that must not be in the file name for the file to be counted. May be left blank to leave files unfiltered.."""
        return len(self._filter_in_dir(file_ext,include,exclude))

    def get_files_in_input_dir(self,file_ext = "",include = "",exclude = ""):
        """Return a list of files in the input directory.

        Parameters:

        file_ext (string, list): the extension of the file type to be returned, or a list of extensions to be returned. May be left blank to include all file extensions.
        
        include (string, list): a character, subtring, or list of characters and/or substrings that must be included in the file name for the file to be returned. May be left blank to leave files unfiltered.
        
        exclude (string, list): a character, subtring, or list of characters and/or substrings that must not be in the file name for the file to be returned. May be left blank to leave files unfiltered."""
        f = self._filter_in_dir(file_ext,include,exclude)
        return f
    
    def get_file_size_in_input_dir(self,file_ext = "",include = "",exclude = "",fmt = "kb"):
        """Return a list of the sizes of files in the input directory.

        Parameters:

        file_ext (string, list): the extension of the file type to determine size of, or a list of extensions to determine size of. May be left blank to include all file extensions.
        
        include (string, list): a character, subtring, or list of characters and/or substrings that must be included in the file name for the file size to be determined. May be left blank to leave files unfiltered.
        
        exclude (string, list): a character, subtring, or list of characters and/or substrings that must not be in the file name for the file size to be determined. May be left blank to leave files unfiltered.
        
        fmt (string): the format of the file sizes to be returned. May be specified as "b" for bytes, "kb" for kilobytes, "mb" for megabytes, or "gb" for gigabytes."""
        def convert_size(size, fmt):
            if fmt == "b": return size
            elif fmt == "kb": return round(size / 1000, 2)
            elif fmt == "mb": return round(size / 1_000_000, 2)
            elif fmt == "gb": return round(size / 1_000_000_000, 2)
            else: raise ValueError(f"Unknown format: {fmt}")

        files = self._filter_in_dir(file_ext, include, exclude)
        return [convert_size(os.stat(os.path.join(self._in_dir, f)).st_size, fmt) for f in files]

    # ============================================================
    # TIMING AND PROCESSING UTILITIES
    # ============================================================

    def initialize_timer(self):
        """Initialize a timer to measure the processing time of a block of code. Returns start time for convenience."""
        from time import perf_counter
        self._start_time = perf_counter()
        return self._start_time

    def return_processing_time(self, decimals = 2):
        """Compute the processing time of a block of code. A timer must have been previously initialized. Returns a tuple of (processing time, units).
        
        Parameters: 
        
        decimals (integer): the number of decimals places the processing time will be rounded to."""
        if self._start_time != None:
            from time import perf_counter
            dur = perf_counter() - self._start_time
            units = "seconds"
            if dur/3600 > 1:
                dur /= 3600
                units = "hours"
            elif dur/60 > 1:
                dur /= 60
                units = "minutes"
            if dur == 1:
                units = units[:-1]
            return (round(dur,decimals), units)
        else:
            print("A timer has not been initialized.")
            return None

    def display_processing_time(self,decimals = 2):
        """Compute and print the processing time of a block of code. A timer must have been previously initialized.
        
        Parameters: 
        
        decimals (integer): the number of decimals places the processing time will be rounded to."""
        if self._start_time != None:
            from time import perf_counter
            dur = perf_counter() - self._start_time
            print(f"Processing completed in {self._format_time(dur,decimals)}.")
        else:
            print("A timer has not been initialized.")

    def initialize_progress_counter(self,total):
        """Initialize a progress counter to provide updates (percent complete) on the progress of iterative processing.
        
        Parameters:
        
        total (integer): the total number of iterations to provide progress updates for."""
        if not isinstance(total, int) or total <= 0:
            raise ValueError("Total iterations must be a positive integer.")
        from time import perf_counter
        self._prog_count = 0
        self._prog_total = total
        self._curr_percent = 0
        self._prog_start_time = perf_counter()

    def update_progress(self, display="percent", timing=None, decimals = 2):
        """Update an initialized progress counter and print the percent completion of an iterative process. Must be called within the iteration loop.
        
        Parameters:

        display (string): the format of the progress output. May be specified as "percent" to print the percent completion, "total" to print the iteration count, or "total_percent" to print both count and percent.
        
        timing (string): optional timing information to include. May be specified as "elapsed", "rate", "eta", "elapsed_rate", "elapsed_eta", "rate_eta", or "elapsed_rate_eta".
        
        decimals (integer): the number of decimals places the processing time will be rounded to."""
        if not isinstance(display, str):
            raise ValueError("'display' must be a string.")
        if display not in {"percent", "total", "total_percent"}:
            raise ValueError("'display' must be one of: 'percent', 'total', 'total_percent'.")

        if timing is not None:
            if not isinstance(timing, str):
                raise ValueError("'timing' must be a string or None.")
            valid_timing = {
                "elapsed", "rate", "eta",
                "elapsed_rate", "elapsed_eta",
                "rate_eta", "elapsed_rate_eta"
            }
            if timing not in valid_timing:
                raise ValueError("'timing' must be one of: " + ", ".join(sorted(valid_timing)) + ".")

        if self._prog_total is None:
            print("Progress counter has not been initialized.")
            return

        from time import perf_counter
        self._prog_count += 1
        percent = int((self._prog_count / self._prog_total) * 100)

        if percent > self._curr_percent:
            self._curr_percent = percent

            if display == "percent":
                output = f"Progress: {percent}%"
            if display == "total":
                output = f"Progress: {self._prog_count} of {self._prog_total}"
            if display == "total_percent":
                output = f"Progress: {self._prog_count} of {self._prog_total} {percent}%"

            if timing is not None:
                raw_elapsed = perf_counter() - self._prog_start_time
                raw_rate = self._prog_count / raw_elapsed if raw_elapsed > 0 else 0
                remaining = self._prog_total - self._prog_count
                raw_eta = remaining / raw_rate if raw_rate > 0 else None

                elapsed = self._format_time(raw_elapsed, decimals)
                eta = self._format_time(raw_eta, decimals) if raw_eta is not None else None

                timing_parts = []
                if "elapsed" in timing:
                    timing_parts.append(f"elapsed: {elapsed}")
                if "rate" in timing:
                    timing_parts.append(f"{round(raw_rate, decimals)} it/s")
                if "eta" in timing:
                    timing_parts.append(f"eta: {eta}" if eta is not None else "eta: --")
                output = output + " | " + " | ".join(timing_parts)

            print(output)

        if self._prog_count >= self._prog_total:
            self._prog_total = None
            self._curr_percent = 0

    def update_progress_bar(self, display="bar", length=50, timing=None, decimals = 2):
        """Update an initialized progress counter and print a progress bar displaying completion of an iterative process. Must be called within the iteration loop.
        
        Parameters:

        display (string): the format of the progress output. May be specified as "bar" to print only the progress bar, "percent" to print the progress bar with percent completion, "total" to print the progress bar with the iteration count, or "total_percent" to print the progress bar with both count and percent.
        
        length (integer): the length of the progress bar in characters. Should be between 5 and 100 characters.
        
        timing (string): optional timing information to include. May be specified as "elapsed", "rate", "eta", "elapsed_rate", "elapsed_eta", "rate_eta", or "elapsed_rate_eta".
        
        decimals (integer): the number of decimals places the processing time will be rounded to."""
        if not isinstance(display, str):
            raise ValueError("'display' must be a string.")
        if display not in {"bar", "percent", "total", "total_percent"}:
            raise ValueError("'display' must be one of: 'bar', 'percent', 'total', 'total_percent'.")
        if not isinstance(length, int):
            raise ValueError("'length' must be an integer.")
        if length < 5 or length > 100:
            raise ValueError("'length' must be between 5 and 100.")

        if timing is not None:
            if not isinstance(timing, str):
                raise ValueError("'timing' must be a string or None.")
            valid_timing = {
                "elapsed", "rate", "eta",
                "elapsed_rate", "elapsed_eta",
                "rate_eta", "elapsed_rate_eta"
            }
            if timing not in valid_timing:
                raise ValueError("'timing' must be one of: " + ", ".join(sorted(valid_timing)) + ".")

        if self._prog_total is None:
            print("Progress counter has not been initialized.")
            return

        from time import perf_counter
        self._prog_count += 1
        percent = int((self._prog_count / self._prog_total) * 100)

        if percent > self._curr_percent:
            self._curr_percent = percent
            num_chevs = int((percent / 100) * length)
            filled = ">" * num_chevs
            empty = "-" * (length - num_chevs)
            update_bar = f"[{filled}{empty}]"

            if display == "bar":
                output = update_bar
            if display == "percent":
                output = f"{update_bar} {percent}%"
            if display == "total":
                output = f"{update_bar} {self._prog_count} of {self._prog_total}"
            if display == "total_percent":
                output = f"{update_bar} {self._prog_count} of {self._prog_total} {percent}%"

            if timing is not None:
                raw_elapsed = perf_counter() - self._prog_start_time
                raw_rate = self._prog_count / raw_elapsed if raw_elapsed > 0 else 0
                remaining = self._prog_total - self._prog_count
                raw_eta = remaining / raw_rate if raw_rate > 0 else None

                elapsed = self._format_time(raw_elapsed, decimals)
                eta = self._format_time(raw_eta, decimals) if raw_eta is not None else None

                timing_parts = []
                if "elapsed" in timing:
                    timing_parts.append(f"elapsed: {elapsed}")
                if "rate" in timing:
                    timing_parts.append(f"{round(raw_rate, decimals)} it/s")
                if "eta" in timing:
                    timing_parts.append(f"eta: {eta}" if eta is not None else "eta: --")
                output = output + " | " + " | ".join(timing_parts)

            print(output, end="\r")

        if self._prog_count >= self._prog_total:
            print()
            self._prog_total = None
            self._curr_percent = 0

    # ============================================================
    # TABLE CREATION AND CONVERSION
    # ============================================================

    def new_table(self,data = None):
        """Create a new Table object to hold tabular data. May be created empty or initialized with an existing Table object ('data').
        
        Parameters:
        
        data (object): the Table object to copy data from. May be left blank to create an empty Table object."""
        t = _Table(data)
        t._boolean_true_values = self._boolean_true_values
        t._boolean_false_values = self._boolean_false_values
        return t
    
    def dict_to_table(self,dict):
        """Attempt to convert a Python dictionary object ('dict') into a Table object and return a new Table.
        
        Parameters:
        
        dict (dictionary): the dictionary to be converted to a Table object. The dictionary keys will be used as column headers, and the dictionary values should be equal-length arrays of values representing column data."""
        # create table object
        t = self.new_table()

        # add columns of data to table
        for h, c in dict.items():
            t.append_column(c,str(h))
        
        # determine data types  
        t.detect_dtypes()

        return t
    
    def table_to_dict(self,table):
        """Attempt to convert a Table object ('table') to a Python dictionary object.
        
        Parameters:
        
        table (object): the Table object to be converted to a Python dictionary. The Table column headers will be used as dictionary keys, and the Table columns will be arrays of values."""
        output = {}
        for h in table.get_headers():
            c = table.get_column(h)
            output[h] = c

        return output

    def load_example_data(self):
        """Load an example dataset for experimentation and demonstration purposes and return a new Table object."""
        package_dir = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.join(package_dir, "Data")
        path = os.path.join(data_dir, "example_data.csv")
        return self.read_csv(path)

    # ============================================================
    # READERS
    # ============================================================

    def peek_raw_text(self,text_file, extension = ".txt", num_lines = 5, return_lines = False):
        """Attempt to read a raw text file ('text_file') and print or return the specified number of lines ('num_lines'). No data type conversion is performed and data is not returned in a Table object. This function streams the file line by line and stops after 'num_lines', so very large files can be inspected safely without loading the entire file into memory.

        Parameters:

        text_file (string): the name of the text file to be read, including the ".txt" extension.
        
        extension (string): the extension of the file to be read.

        num_lines (integer): the number of lines to be printed or returned.
        
        return_lines (Boolean): flag to indicate if num_lines will be printed to the terminal (False) or returned as a list (True)."""
        # open file
        file = self._extension_handler(text_file,extension)
        if self._in_dir and os.path.dirname(file) == "":
            file = os.path.join(self._in_dir, file)

        lines = []
        with open(file, "r", encoding="utf-8") as f:
            for i, line in enumerate(f):
                if not return_lines:
                    print(line.rstrip("\n"))
                else:
                    lines.append(line.rstrip("\n"))
                if i+1 == num_lines:
                    if return_lines:
                        return lines
                    break

    def read_raw_text(self,text_file, extension = ".txt"):
        """Attempt to read a raw text file ('text_file') and return the contents as a nested list. No data type conversion is performed and data is not returned in a Table object.

        Parameters:

        text_file (string): the name of the text file to be read, including the ".txt" extension.
        
        extension (string): the extension of the file to be read."""
        # open file
        file = self._extension_handler(text_file,extension)
        if self._in_dir and os.path.dirname(file) == "":
            file = os.path.join(self._in_dir, file)
        
        with open(file, "r", encoding="utf-8") as f:
            return [line.rstrip("\n") for line in f.readlines()]

    def read_text(self,text_file,delimiter=",",skiprows = 0,headers = True):
        """Attempt to read a delimited text file ('text_file') with a ".txt" extension and store the data in a new Table object.

        Parameters:
        
        text_file (string): the name of the text file to be read, including the ".txt" extension.
        
        delimiter (string): the character that separates values in the text file. A comma is the default delimiter.
        
        skiprows (integer): the number of rows to skip at the start of the file, if necessary.
        
        headers (Boolean): flag to indicate if the columns have headers. If columns do not have headers, generic headers will be automatically generated."""
        file = self._extension_handler(text_file,".txt")
        if self._in_dir and os.path.dirname(file) == "":
            file = os.path.join(self._in_dir, file)
        
        rows = []
        with open(file, "r", encoding="utf-8") as f:
            for row in self._read_logical_rows(f):
                # fast path: no quotes, use split (very fast)
                if '"' not in row:
                    rows.append(row.split(delimiter))
                else:
                    # slow path: quoted fields, use robust splitter
                    rows.append(self._split_quoted_line(row, delimiter=delimiter))

        # create table object
        t = self.new_table()
        
        # skiprows, if necessary
        if skiprows > 0:
            rows = rows[skiprows:]
        
        # headers
        if headers:
            col_headers = rows[0]
            data = [self._row_dtype(row) for row in rows[1:]]
            for h in range(len(col_headers)):#remove extra quotes from headers
                if col_headers[h].endswith('"') and col_headers[h].startswith('"'): 
                    col_headers[h] = col_headers[h][1:-1]
            t._load_on_read(data)
            t._headers = col_headers
        else:
            data = [self._row_dtype(row) for row in rows]
            t._load_on_read(data)
            t.autogenerate_headers()
        
        # determine data types
        t.detect_dtypes()
        
        return t

    def read_other_text(self,text_file,extension,delimiter=",",skiprows = 0,headers = True):
        """Attempt to read a delimited text file ('text_file') with an extension other than ".txt" and store the data in a new Table object.

        Parameters:
        
        text_file (string): the name of the text file to be read.

        extension (string): the extension of the file to be read.
        
        delimiter (string): the character that separates values in the text file. A comma is the default delimiter.
        
        skiprows (integer): the number of rows to skip at the start of the file, if necessary.
        
        headers (Boolean): flag to indicate if the columns have headers. If columns do not have headers, generic headers will be automatically generated."""
        # open file
        file = self._extension_handler(text_file,extension)
        if self._in_dir and os.path.dirname(file) == "":
            file = os.path.join(self._in_dir, file)

        rows = []
        with open(file, "r", encoding="utf-8") as f:
            for row in self._read_logical_rows(f):
                # fast path: no quotes, use split (very fast)
                if '"' not in row:
                    rows.append(row.split(delimiter))
                else:
                    # slow path: quoted fields, use robust splitter
                    rows.append(self._split_quoted_line(row, delimiter=delimiter))

        # create table object
        t = self.new_table()

        # skiprows, if necessary
        if skiprows > 0:
            rows = rows[skiprows:]
        
        # headers
        if headers:
            col_headers = rows[0]
            data = [self._row_dtype(row) for row in rows[1:]]
            for h in range(len(col_headers)):#remove extra quotes from headers
                if col_headers[h].endswith('"') and col_headers[h].startswith('"'): 
                    col_headers[h] = col_headers[h][1:-1]
            t._load_on_read(data)
            t._headers = col_headers
        else:
            data = [self._row_dtype(row) for row in rows]
            t._load_on_read(data)
            t.autogenerate_headers()
        
        # determine data types
        t.detect_dtypes()
        
        return t
    
    def read_csv(self,csv_file,skiprows = 0,headers = True):
        """Attempt to read a delimited CSV file ('csv_file') and store the data in a new Table object.

        Parameters:
        
        csv_file (string): the name of the CSV file to be read, including the ".csv" extension.
        
        skiprows (integer): the number of rows to skip at the start of the file, if necessary.
        
        headers (Boolean): flag to indicate if the columns have headers. If columns do not have headers, generic headers will be automatically generated."""
        # open file
        file = self._extension_handler(csv_file,".csv")
        if self._in_dir and os.path.dirname(file) == "":
            file = os.path.join(self._in_dir, file)

        rows = []
        with open(file, "r", encoding="utf-8") as f:
            for row in self._read_logical_rows(f):
                # fast path: no quotes, use split (very fast)
                if '"' not in row:
                    rows.append(row.split(","))
                else:
                    # slow path: quoted fields, use robust splitter
                    rows.append(self._split_quoted_line(row, delimiter=","))

        # create table object
        t = self.new_table()

        # skiprows, if necessary
        if skiprows > 0:
            rows = rows[skiprows:]
        
        # headers
        if headers:
            col_headers = rows[0]
            data = [self._row_dtype(row) for row in rows[1:]]
            for h in range(len(col_headers)):#remove extra quotes from headers
                if col_headers[h].endswith('"') and col_headers[h].startswith('"'): 
                    col_headers[h] = col_headers[h][1:-1]
            t._load_on_read(data)
            t._headers = col_headers
        else:
            data = [self._row_dtype(row) for row in rows]
            t._load_on_read(data)
            t.autogenerate_headers()
        
        # determine data types
        t.detect_dtypes()

        return t

    def read_tsv(self, tsv_file, skiprows=0, headers=True):
        """Attempt to read a tab-separated values file ('tsv_file') and store the data in a new Table object.

        Parameters:

        tsv_file (string): the name of the TSV file to be read, including the ".tsv" extension.

        skiprows (integer): the number of rows to skip at the start of the file, if necessary.

        headers (Boolean): flag to indicate if the columns have headers. If columns do not have headers, generic headers will be automatically generated."""
        # Ensure extension and working directory
        file = self._extension_handler(tsv_file,".tsv")
        if self._in_dir and os.path.dirname(file) == "":
            file = os.path.join(self._in_dir, file)

        rows = []
        with open(file, "r", encoding="utf-8") as f:
            for row in self._read_logical_rows(f):
                # fast path: no quotes, use split (very fast)
                if '"' not in row:
                    rows.append(row.split("\t"))
                else:
                    # slow path: quoted fields, use robust splitter
                    rows.append(self._split_quoted_line(row, delimiter="\t"))

        # Create table object
        t = self.new_table()

        # Skiprows
        if skiprows > 0:
            rows = rows[skiprows:]

        # Handle headers
        if headers and rows:
            col_headers = rows[0]
            data = [self._row_dtype(row) for row in rows[1:]]
            t._load_on_read(data)
            t._headers = col_headers
        else:
            data = [self._row_dtype(row) for row in rows]
            t._load_on_read(data)
            t.autogenerate_headers()

        # Detect column data types
        t.detect_dtypes()

        return t

    def read_psv(self, psv_file, skiprows=0, headers=True):
        """Attempt to read a pipe-separated values file ('psv_file') and store the data in a new Table object.

        Parameters:

        psv_file (string): the name of the PSV file to be read, including the ".psv" extension.

        skiprows (integer): the number of rows to skip at the start of the file, if necessary.

        headers (Boolean): flag to indicate if the columns have headers. If columns do not have headers, generic headers will be automatically generated."""
        # Ensure extension and working directory
        file = self._extension_handler(psv_file,".psv")
        if self._in_dir and os.path.dirname(file) == "":
            file = os.path.join(self._in_dir, file)

        # Read raw lines
        rows = []
        with open(file, "r", encoding="utf-8") as f:
            for row in self._read_logical_rows(f):
                # fast path: no quotes, use split (very fast)
                if '"' not in row:
                    rows.append(row.split("|"))
                else:
                    # slow path: quoted fields, use robust splitter
                    rows.append(self._split_quoted_line(row, delimiter="|"))

        # Create table object
        t = self.new_table()

        # Skiprows
        if skiprows > 0:
            rows = rows[skiprows:]

        # Handle headers
        if headers and rows:
            col_headers = rows[0]
            data = [self._row_dtype(row) for row in rows[1:]]
            t._load_on_read(data)
            t._headers = col_headers
        else:
            data = [self._row_dtype(row) for row in rows]
            t._load_on_read(data)
            t.autogenerate_headers()

        # Detect column data types
        t.detect_dtypes()

        return t

    def read_fixedwidth(self, fw_file, extension = ".txt", widths=None, skiprows=0, headers=True, overflow=None):
        """Attempt to read a fixed-width text file ('txt_file') and store the data in a new Table object.

        Parameters:

        fw_file (string): the name of the fixed-width file to be read, including the extension.

        extension (string): the extension of the file to be read.

        widths (list): a list of character widths of each column. If None, column boundaries will be auto-detected based on whitespace patterns.

        skiprows (integer): the number of rows to skip at the start of the file, if necessary.

        headers (Boolean): flag to indicate if the columns have headers. If columns do not have headers, generic headers will be automatically generated.

        overflow (string): how column overflow should be managed if a column value exceeds its specified width. If None, columns are silently truncationed. If "auto-expand", columns are expanded as need to fit values. If "strict", an error is raised."""

        # Ensure extension and working directory
        file = self._extension_handler(fw_file,extension)
        if self._in_dir and os.path.dirname(file) == "":
            file = os.path.join(self._in_dir, file)

        # Read raw lines
        with open(file, "r", encoding="utf-8") as f:
            lines = [line.rstrip("\n") for line in f.readlines()]

        # Skip initial non-tabular lines if requested
        if skiprows > 0:
            lines = lines[skiprows:]

        # Auto-detect widths if not provided
        if widths is None:
            sample = next((ln for ln in lines if ln.strip()), "")
            boundaries = []
            in_field = False

            for i, ch in enumerate(sample):
                if not in_field and ch != " ":
                    boundaries.append(i)
                    in_field = True
                elif in_field and ch == " ":
                    in_field = False

            widths = []
            for i in range(len(boundaries)):
                start = boundaries[i]
                end = boundaries[i + 1] if i + 1 < len(boundaries) else len(sample)
                widths.append(end - start)

        # Parse rows using widths + overflow behavior
        rows = []
        for line in lines:
            if not line.strip():
                continue

            row = []
            pos = 0

            for i, w in enumerate(widths):
                # Basic slice
                cell = line[pos:pos + w]

                # STRICT MODE
                if overflow == "strict":
                    if len(cell.rstrip()) > w:
                        raise ValueError(
                            f"Fixed-width overflow detected in column {i+1}: "
                            f"'{cell.rstrip()}' exceeds width {w}"
                        )
                    row.append(cell.strip())
                    pos += w
                    continue

                # AUTO-EXPAND MODE
                if overflow == "auto-expand":
                    # Look ahead to next field start
                    if i < len(widths) - 1:
                        next_start = pos + w

                        # Overflow detected: next field does NOT start with space
                        if next_start < len(line) and line[next_start] != " ":
                            j = next_start
                            # Extend until whitespace or end of line
                            while j < len(line) and line[j] != " ":
                                cell += line[j]
                                j += 1
                            pos = j  # move to new boundary
                            row.append(cell.strip())
                            continue

                    # No overflow, normal slice
                    row.append(cell.strip())
                    pos += w
                    continue

                # DEFAULT MODE (overflow=None)
                row.append(cell.strip())
                pos += w

            rows.append(row)

        # Create table object
        t = self.new_table()

        if not rows:
            t._load_on_read([])
            t.autogenerate_headers()
            return t

        # Handle headers
        if headers:
            col_headers = rows[0]
            data_rows = rows[1:]
            data = [self._row_dtype(row) for row in data_rows]
            t._load_on_read(data)
            t._headers = col_headers
        else:
            data = [self._row_dtype(row) for row in rows]
            t._load_on_read(data)
            t.autogenerate_headers()

        # Detect column data types
        t.detect_dtypes()

        return t

    def read_markdown(self, md_file, headers=True, skiprows=0, table_index=0):
        """Attempt to read a table in a Markdown document ('md_file') and store the data in a new Table object.

        Parameters:

        md_file (string): the name of the Markdown file to be read, including the ".md" extension.

        headers (Boolean): flag to indicate if the columns have headers. If columns do not have headers, generic headers will be automatically generated.

        skiprows (integer): the number of rows to skip at the start of the file, if necessary.

        table_index (integer): which table to read. Markdown documents may contain multiple tables. Index 0 refers to the first table found."""
        # Ensure extension and working directory
        file = self._extension_handler(md_file,".md",valid_exts=[".md", ".markdown"])
        if self._in_dir and os.path.dirname(file) == "":
            file = os.path.join(self._in_dir, file)

        # Read raw lines
        with open(file, "r", encoding="utf-8") as f:
            lines = [line.rstrip("\n") for line in f.readlines()]

        # Skip initial non-tabular lines if requested
        if skiprows > 0:
            lines = lines[skiprows:]

        # Markdown table detection:
        # A table is defined as:
        #   1. A header row starting with '|'
        #   2. A separator row containing '-' and '|'
        #   3. One or more data rows starting with '|'
        tables = []
        current_table = []
        in_table = False
        saw_separator = False

        for line in lines:
            stripped = line.strip()

            # Detect header row
            if stripped.startswith("|") and "|" in stripped:
                # If we were not in a table, start a new one
                if not in_table:
                    in_table = True
                    saw_separator = False
                    current_table = []

                # Add row if it's not a separator
                if not (set(stripped.replace("|", "").replace(":", "").strip()) <= {"-"}):
                    # Normal row
                    cells = [cell.strip() for cell in stripped.strip("|").split("|")]
                    current_table.append(cells)
                else:
                    # Separator row
                    saw_separator = True

            else:
                # If we were in a table and hit a non-table line, close the table
                if in_table:
                    if saw_separator and current_table:
                        tables.append(current_table)
                    in_table = False
                    current_table = []
                    saw_separator = False

        # Edge case: file ends while still in a table
        if in_table and saw_separator and current_table:
            tables.append(current_table)

        # Select the requested table
        rows = tables[table_index] if table_index < len(tables) else []

        # Create table object
        t = self.new_table()

        if not rows:
            # No table found, return empty table with autogenerated headers
            t._load_on_read([])
            t.autogenerate_headers()
            return t

        # Handle headers
        if headers:
            col_headers = rows[0]
            data_rows = rows[1:]
            data = [self._row_dtype(row) for row in data_rows]
            t._load_on_read(data)
            t._headers = col_headers
        else:
            data = [self._row_dtype(row) for row in rows]
            t._load_on_read(data)
            t.autogenerate_headers()

        # Detect column data types
        t.detect_dtypes()

        return t

    def read_latex(self, tex_file, headers=True, skiprows=0, table_index=0):
        """Attempt to read a table in a LaTeX document ('tex_file') and store the data in a new Table object.

        Parameters:

        tex_file (string): the name of the LaTeX file to be read, including the ".tex" extension.

        headers (Boolean): flag to indicate if the columns have headers. If columns do not have headers, generic headers will be automatically generated.

        skiprows (integer): the number of rows to skip at the start of the file, if necessary.

        table_index (integer): which table to read. LaTeX documents may contain multiple tables. Index 0 refers to the first table found."""
        # Ensure extension and working directory
        file = self._extension_handler(tex_file,".tex")
        if self._in_dir and os.path.dirname(file) == "":
            file = os.path.join(self._in_dir, file)

        # Read raw lines
        with open(file, "r", encoding="utf-8") as f:
            lines = [line.rstrip("\n") for line in f.readlines()]

        # Skip initial non-tabular lines if requested
        if skiprows > 0:
            lines = lines[skiprows:]

        # Detect LaTeX tables:
        #   \begin{tabular}{...}
        #   \begin{longtable}{...}
        #   ... rows ...
        #   \end{tabular}
        #   \end{longtable}
        tables = []
        current_table = []
        in_table = False

        for line in lines:
            stripped = line.strip()

            # Remove comments
            if "%" in stripped:
                stripped = stripped.split("%", 1)[0].strip()

            # Detect table start
            if stripped.startswith(r"\begin{tabular") or stripped.startswith(r"\begin{longtable"):
                in_table = True
                current_table = []
                continue

            # Detect table end
            if stripped.startswith(r"\end{tabular") or stripped.startswith(r"\end{longtable"):
                if current_table:
                    tables.append(current_table)
                in_table = False
                current_table = []
                continue

            if in_table:
                # Skip formatting commands
                if stripped in ("\\hline", r"\hline", ""):
                    continue

                # Detect row endings
                if "\\" in stripped:
                    row_text = stripped.split("\\")[0].strip()
                    cells = [cell.strip() for cell in row_text.split("&")]
                    current_table.append(cells)

        # Select the requested table
        rows = tables[table_index] if table_index < len(tables) else []

        # Create table object
        t = self.new_table()

        if not rows:
            # No table found, return empty table with autogenerated headers
            t._load_on_read([])
            t.autogenerate_headers()
            return t

        # Handle headers
        if headers:
            col_headers = rows[0]
            data_rows = rows[1:]
            data = [self._row_dtype(row) for row in data_rows]
            t._load_on_read(data)
            t._headers = col_headers
        else:
            data = [self._row_dtype(row) for row in rows]
            t._load_on_read(data)
            t.autogenerate_headers()

        # Detect column data types
        t.detect_dtypes()

        return t

    def read_yaml(self, yaml_file):
        """Attempt to read a YAML document ('yaml_file') containing a list of dictionaries and store the data in a new Table object.

        The YAML file must contain a top-level list, where each element is a dictionary representing a row. Dictionary keys become column headers.

        Parameters:

        yaml_file (string): the name of the YAML file to be read, including the ".yaml" or ".yml" extension."""
        # Ensure extension and working directory
        file = self._extension_handler(yaml_file, ".yaml", valid_exts=[".yaml", ".yml"])
        if self._in_dir and os.path.dirname(file) == "":
            file = os.path.join(self._in_dir, file)

        # Read raw lines
        with open(file, "r", encoding="utf-8") as f:
            lines = [line.rstrip("\n") for line in f.readlines()]

        # Minimal YAML parser for list-of-dicts
        rows = []
        current_dict = None

        for line in lines:
            stripped = line.strip()

            # Skip empty lines and comments
            if not stripped or stripped.startswith("#"):
                continue

            # Start of a new dictionary
            if stripped.startswith("- "):
                if current_dict is not None:
                    rows.append(current_dict)
                current_dict = {}
                stripped = stripped[2:].strip()  # remove "- "

                # Handle inline key-value after "- "
                if ":" in stripped:
                    key, value = stripped.split(":", 1)
                    current_dict[key.strip()] = value.strip()

            # Key-value inside a dictionary
            elif ":" in stripped and current_dict is not None:
                key, value = stripped.split(":", 1)
                current_dict[key.strip()] = value.strip()

        # Append last dictionary
        if current_dict is not None:
            rows.append(current_dict)

        # Create table object
        t = self.new_table()

        if not rows:
            # No data, empty table
            t._load_on_read([])
            t.autogenerate_headers()
            return t

        # Determine all headers (union of keys)
        headers = sorted({key for row in rows for key in row.keys()})

        # Build row lists in header order
        data = []
        for row in rows:
            ordered = [row.get(h, None) for h in headers]
            data.append(self._row_dtype(ordered))

        # Load into table
        t._load_on_read(data)
        t._headers = headers

        # Detect column data types
        t.detect_dtypes()

        return t

    def read_ini(self, ini_file, headers=True, skiprows=0):
        """Attempt to read an INI-style configuration file ('ini_file') and store the data in a new Table object. Each section becomes a row, and keys within sections become columns. Missing keys are filled with None.

        Parameters:

        ini_file (string): the name of the INI file to be read. If no file extension is provided, '.ini' will be appended automatically.

        headers (Boolean): flag to indicate if the keys should be used as column headers. If False, generic headers will be automatically generated.

        skiprows (integer): the number of rows to skip at the start of the file, if necessary."""
        # Ensure extension and working directory
        file = self._extension_handler(ini_file, ".ini", valid_exts=[".ini", ".cfg", ".conf"])
        if self._in_dir and os.path.dirname(file) == "":
            file = os.path.join(self._in_dir, file)

        # Read raw lines
        with open(file, "r", encoding="utf-8") as f:
            lines = [line.rstrip("\n") for line in f.readlines()]

        # Skip initial non-tabular lines if requested
        if skiprows > 0:
            lines = lines[skiprows:]

        rows = []
        current_section = None
        current_dict = {}

        for line in lines:
            stripped = line.strip()

            # Skip empty lines and comments
            if not stripped or stripped.startswith(";") or stripped.startswith("#"):
                continue

            # Section header: [SectionName]
            if stripped.startswith("[") and stripped.endswith("]"):
                # Save previous section
                if current_section is not None:
                    rows.append(current_dict)
                current_section = stripped[1:-1].strip()
                current_dict = {}
                continue

            # Key-value pair: key = value
            if "=" in stripped and current_section is not None:
                key, value = stripped.split("=", 1)
                key = key.strip()
                value = value.strip()
                current_dict[key] = value
                continue

        # Append last section
        if current_section is not None:
            rows.append(current_dict)

        # Create table object
        t = self.new_table()

        if not rows:
            t._load_on_read([])
            t.autogenerate_headers()
            return t

        # Determine all headers (union of keys)
        all_headers = sorted({key for row in rows for key in row.keys()})

        # Build row lists in header order
        data = []
        for row in rows:
            ordered = [row.get(h, None) for h in all_headers]
            data.append(self._row_dtype(ordered))

        # Load into table
        t._load_on_read(data)

        if headers:
            t._headers = all_headers
        else:
            t.autogenerate_headers()

        # Detect column data types
        t.detect_dtypes()

        return t

    def read_xml(self, xml_file, headers=True, skiprows=0, table_index=0):
        """Attempt to read a table in an XML document ('xml_file') and store the data in a new Table object.

        Parameters:

        xml_file (string): the name of the XML file to be read, including the ".xml" extension.

        headers (Boolean): flag to indicate if the columns have headers. If columns do not have headers, generic headers will be automatically generated.

        skiprows (integer): the number of rows to skip at the start of the file, if necessary.

        table_index (integer): which table to read. XML documents may contain multiple <table> elements. Index 0 refers to the first table found."""
        # Ensure extension and working directory
        file = self._extension_handler(xml_file, ".xml")
        if self._in_dir and os.path.dirname(file) == "":
            file = os.path.join(self._in_dir, file)

        # Read raw lines
        with open(file, "r", encoding="utf-8") as f:
            lines = [line.rstrip("\n") for line in f.readlines()]

        # Skip initial non-tabular lines if requested
        if skiprows > 0:
            lines = lines[skiprows:]

        # Parse XML manually (zero dependencies)
        tables = []
        current_table = []
        in_table = False
        in_row = False
        current_row = {}

        for line in lines:
            stripped = line.strip()

            # Skip comments
            if stripped.startswith("<!--"):
                continue

            # Detect <table>
            if stripped.startswith("<table"):
                in_table = True
                current_table = []
                continue

            # Detect </table>
            if stripped.startswith("</table>"):
                if current_table:
                    tables.append(current_table)
                in_table = False
                current_table = []
                continue

            if not in_table:
                continue

            # Detect <row>
            if stripped.startswith("<row"):
                in_row = True
                current_row = {}
                continue

            # Detect </row>
            if stripped.startswith("</row>"):
                if current_row:
                    current_table.append(current_row)
                in_row = False
                current_row = {}
                continue

            # Inside a row: extract <tag>value</tag>
            if in_row and "<" in stripped and "</" in stripped:
                # Extract all tags on the line
                pos = 0
                while True:
                    start_tag = stripped.find("<", pos)
                    if start_tag == -1:
                        break
                    end_tag = stripped.find(">", start_tag)
                    if end_tag == -1:
                        break

                    tag = stripped[start_tag + 1:end_tag].strip()
                    close_tag = f"</{tag}>"
                    close_pos = stripped.find(close_tag, end_tag)
                    if close_pos == -1:
                        break

                    value = stripped[end_tag + 1:close_pos].strip()
                    current_row[tag] = value

                    pos = close_pos + len(close_tag)

        # Select the requested table
        rows = tables[table_index] if table_index < len(tables) else []

        # Create table object
        t = self.new_table()

        if not rows:
            t._load_on_read([])
            t.autogenerate_headers()
            return t

        # Determine headers (union of all tags)
        all_headers = sorted({key for row in rows for key in row.keys()})

        # Build row lists
        data = []
        for row in rows:
            ordered = [row.get(h, None) for h in all_headers]
            data.append(self._row_dtype(ordered))

        # Load into table
        t._load_on_read(data)

        if headers:
            # Use tag names as headers
            t._headers = all_headers
        else:
            # Treat tag names as data; autogenerate headers
            t.autogenerate_headers()

        # Detect column data types
        t.detect_dtypes()

        return t

    def read_json(self,json_file):
        """Attempt to read a JSON file ('json_file') and store the data in a new Table object. It is expected that the JSON file has a single object. The object's keys will be used as column headers, and the object's values should be equal-length arrays representing the column data.
        
        Parameters:
        
        json_file (string): the name of the JSON file to be read, including the ".json" extension."""
        import json
        file = self._extension_handler(json_file, ".json")
        if self._in_dir and os.path.dirname(file) == "":
            file = os.path.join(self._in_dir, file)

        # parse JSON file json module
        with open(file, "r", encoding="utf-8") as f:
            jdict = json.load(f)

        # create table object
        t = self.new_table()

        # add columns of data to table
        for h, c in jdict.items():
            t.append_column(c,str(h))
        
        # determine data types  
        t.detect_dtypes()

        return t

    def read_html(self,html_file,headers = True, table_index = 0):
        """Attempt to read a table in an HTML document ('html_file') and store the data in a new Table object.

        Parameters:
        
        html_file (string): the name of the HTML file to be read, including the ".html" extension.
        
        headers (Boolean): flag to indicate if the columns have headers. If columns do not have headers, generic headers will be automatically generated.
        
        table_index (integer): which table to read. In an HTML document with one table, index must be 0 (first table). In a document with multiple tables, other index values may be used to read other tables."""
        # check extension and working directory, then open the file
        file = self._extension_handler(html_file, ".html", valid_exts=[".html", ".htm"])
        if self._in_dir and os.path.dirname(file) == "":
            file = os.path.join(self._in_dir, file)

        with open(file, "r", encoding="utf-8") as f:
            lines = f.readlines()

        # create table object
        t = self.new_table()
        tables = []
        current_table = []
        current_row = []
        in_table = False

        for line in lines:
            line = line.strip()
            if "<table" in line:
                in_table = True
                current_table = []
            elif "</table>" in line and in_table:
                if current_row:
                    current_table.append(current_row)
                    current_row = []
                tables.append(current_table)
                in_table = False
            elif in_table:
                # extract ALL <th> cells in the line
                while "<th>" in line and "</th>" in line:
                    start = line.find("<th>") + 4
                    end = line.find("</th>", start)
                    cell = line[start:end].strip()
                    current_row.append(cell)
                    line = line[end+5:]  # move past this cell

                # extract ALL <td> cells in the line
                while "<td>" in line and "</td>" in line:
                    start = line.find("<td>") + 4
                    end = line.find("</td>", start)
                    cell = line[start:end].strip()
                    current_row.append(cell)
                    line = line[end+5:]  # move past this cell

                # row boundary
                if "</tr>" in line and current_row:
                    current_table.append(current_row)
                    current_row = []

        rows = tables[table_index] if table_index < len(tables) else []

        if headers and rows:
            col_headers = rows[0]
            data = [self._row_dtype(row) for row in rows[1:]]
            for i, h in enumerate(col_headers):
                if h.startswith('"') and h.endswith('"'):
                    col_headers[i] = h[1:-1]
            t._load_on_read(data)
            t._headers = col_headers
        else:
            data = [self._row_dtype(row) for row in rows]
            t._load_on_read(data)
            t.autogenerate_headers()

        t.detect_dtypes()
        return t

    def read_r_dput(self, r_file, headers=True, skiprows=0):
        """Attempt to read an R dput() representation of a data frame ('r_file') and store the data in a new Table object.

        The R object must be of the form:
            structure(list(
                col1 = c(...),
                col2 = c(...),
                ...
            ))

        Parameters:

        r_file (string): the name of the R dput file to be read. If no file extension is provided, '.r' will be appended automatically.

        headers (Boolean): flag to indicate if the column names should be used as headers. If False, generic headers will be automatically generated.

        skiprows (integer): the number of rows to skip at the start of the file, if necessary."""
        # Ensure extension and working directory
        file = self._extension_handler(r_file,  ".r", valid_exts=[".r", ".dput"])
        if self._in_dir and os.path.dirname(file) == "":
            file = os.path.join(self._in_dir, file)

        # Read raw lines
        with open(file, "r", encoding="utf-8") as f:
            lines = [line.rstrip("\n") for line in f.readlines()]

        # Skip initial non-tabular lines if requested
        if skiprows > 0:
            lines = lines[skiprows:]

        # Join into a single string for easier parsing
        text = " ".join(lines)

        # Extract the list(...) content
        if "list(" not in text:
            raise ValueError("File does not contain an R list(...) structure.")

        start = text.find("list(") + 5
        end = text.rfind(")")
        content = text[start:end].strip()

        # Split top-level list entries: colname = c(...)
        columns = {}
        depth = 0
        token = ""
        entries = []

        for ch in content:
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1

            if ch == "," and depth == 0:
                entries.append(token.strip())
                token = ""
            else:
                token += ch

        if token.strip():
            entries.append(token.strip())

        # Parse each entry: name = c(...)
        for entry in entries:
            if "=" not in entry:
                continue

            name, vec = entry.split("=", 1)
            name = name.strip()

            # Extract c(...)
            vec = vec.strip()
            if vec.startswith("c(") and vec.endswith(")"):
                vec = vec[2:-1].strip()

            # Split vector elements
            values = []
            current = ""
            in_quotes = False

            for ch in vec:
                if ch == '"' and not in_quotes:
                    in_quotes = True
                    current += ch
                elif ch == '"' and in_quotes:
                    in_quotes = False
                    current += ch
                elif ch == "," and not in_quotes:
                    values.append(current.strip())
                    current = ""
                else:
                    current += ch

            if current.strip():
                values.append(current.strip())

            # Clean values: remove quotes, convert NA, numeric, logical
            cleaned = []
            for v in values:
                v = v.strip()

                if v == "NA":
                    cleaned.append(None)
                elif v.startswith('"') and v.endswith('"'):
                    cleaned.append(v[1:-1])
                elif v.lower() in ("true", "false"):
                    cleaned.append(v.lower() == "true")
                else:
                    # Try numeric conversion
                    try:
                        cleaned.append(float(v) if "." in v else int(v))
                    except:
                        cleaned.append(v)

            columns[name] = cleaned

        # Determine number of rows
        lengths = {len(v) for v in columns.values()}
        if len(lengths) != 1:
            raise ValueError("Column vectors have inconsistent lengths.")

        nrows = lengths.pop()

        # Build row-wise data
        data = []
        col_names = sorted(columns.keys())

        for i in range(nrows):
            row = [columns[c][i] for c in col_names]
            data.append(self._row_dtype(row))

        # Create table object
        t = self.new_table()
        t._load_on_read(data)

        if headers:
            t._headers = col_names
        else:
            t.autogenerate_headers()

        t.detect_dtypes()
        return t

    def read_sql(self,sql_file,table_name):
        """Attempt to read a table ('table_name') from the specified SQL database ('sql_file') and store the data in a new Table object.

        Parameters:
        
        sql_file (string): the name of the SQL database file to be read, including one of ".db", ".db2", ".db3", ".sqlite", ".sqlite2", or ".sqlite3" as an extension. If no extension is provided, the function will attempt to open a file with the ".sqlite" extension.
        
        table_name (string): the name of the table to be read from the specified database file. """
        # open the file
        import sqlite3
        sql_file = self._extension_handler(sql_file,target_ext=".sqlite", valid_exts=[".db", ".db2", ".db3", ".sqlite", ".sqlite2", ".sqlite3"])
        if self._in_dir and os.path.dirname(sql_file) == "":
            sql_file = os.path.join(self._in_dir, sql_file)

        # create table object
        t = self.new_table()

        # open connection safely
        with sqlite3.connect(sql_file) as conn:
            conn.row_factory = sqlite3.Row
            crs = conn.cursor()
            crs.execute(f"SELECT * FROM {table_name}")
            rows = crs.fetchall()

            if rows:
                headers = list(rows[0].keys())
                data = [self._row_dtype(list(row)) for row in rows]

                # clean headers
                for i, h in enumerate(headers):
                    if h.startswith('"') and h.endswith('"'):
                        headers[i] = h[1:-1]

                t._load_on_read(data)
                t._headers = headers

        # determine data types
        t.detect_dtypes()

        return t

    def read_dbf(self,dbf_file):
        """Attempt to read a Dbase or Xbase database file ('dbf_file') with a ".dbf" extension and store the data in a Table object. Null values in the file will be read as empty strings ('').

        Parameters:
            
        dbf_file (string): the name of the Dbase or Xbase file to be read, including the ".dbf" extension."""
        # The core components of this Dbase/Xbase reader (unpacking and decoding the binary data) are based on the works of:
        # Raymond Hettinger (https://code.activestate.com/recipes/362715-dbf/)
        # Tomas Nordin (https://code.activestate.com/recipes/580696-dbf-reader-and-writer-selective-fields-and-nullrep/)
        # The first two rows contain the field names and specifications (type, size, decimal places), respectively
        # The following rows contain the data.

        # open the file
        import struct, decimal
        file = self._extension_handler(dbf_file, ".dbf")
        if self._in_dir and os.path.dirname(file) == "":
            file = os.path.join(self._in_dir, file)

        # create table object
        t = self.new_table()

        with open(file, "rb") as f:
            # read header
            num_rows, len_header = struct.unpack("<xxxxLH22x", f.read(32))
            num_cols = (len_header - 33) // 32

            headers = []
            for _ in range(num_cols):
                name, typ, size, dec = struct.unpack("<11sc4xBB14x", f.read(32))
                name = name.decode("ascii", errors="ignore").replace("\0", "")
                typ = typ.decode("ascii", errors="ignore")
                headers.append((name, typ, size, dec))

            output = [[col[0] for col in headers]]

            # skip header terminator
            f.read(1)

            # add deletion flag
            headers.insert(0, ("DeletionFlag", "C", 1, 0))
            fmt = "".join(["%ds" % f_i[2] for f_i in headers])
            fmt_size = struct.calcsize(fmt)

            for _ in range(num_rows):
                recordb = struct.unpack(fmt, f.read(fmt_size))
                record = [field.decode("ascii", errors="ignore") for field in recordb]

                # skip deleted records
                if not record or record[0] != " ":
                    continue

                result = []
                for (name, typ, size, deci), value in zip(headers, record):
                    if name == "DeletionFlag":
                        continue
                    if typ == "N":  # number
                        value = value.replace("\0", "").strip()
                        if value == "" or "*" in value:
                            value = ""
                        elif deci:
                            value = decimal.Decimal(value)
                        else:
                            value = int(value)
                    elif typ == "D":  # date
                        if len(value) >= 8:
                            y, m, d = int(value[:4]), int(value[4:6]), int(value[6:8])
                            value = datetime.date(y, m, d)
                    elif typ == "L":  # logical
                        value = ("T" if value in "YyTt" else
                                "F" if value in "NnFf" else "?")
                    elif typ == "F":  # float
                        try:
                            value = float(value)
                        except ValueError:
                            pass
                    result.append(value)
                output.append([str(v).rstrip() for v in result])

        rows = output

        # headers
        col_headers = rows[0]
        data = [self._row_dtype(row) for row in rows[1:]]
        for i, h in enumerate(col_headers):
            if h.startswith('"') and h.endswith('"'):
                col_headers[i] = h[1:-1]
        t._load_on_read(data)
        t._headers = col_headers

        # determine data types
        t.detect_dtypes()

        return t

    def read_xlsx(self, xlsx_file, sheet_name, headers = True):
        """Attempt to read a sheet ('sheet_name') from the specified Excel workbook ('xlsx_file') and store the data in a new Table object.

        Parameters:
        
        xlsx_file (string): the name of the Excel workbook file to be read, including the ".xlsx" extension.
        
        sheet_name (string): the name of the sheet to be read from the specified workbook file.
        
        headers (Boolean): flag to indicate if the columns have headers. If columns do not have headers, generic headers will be automatically generated."""
        import zipfile
        import xml.etree.ElementTree as ET

        NS = {
            "main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
            "rel": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
            "pkgrel": "http://schemas.openxmlformats.org/package/2006/relationships",
        }

        def _col_letters_to_index(col_letters: str) -> int:
            idx = 0
            for ch in col_letters:
                idx = idx * 26 + (ord(ch) - ord("A") + 1)
            return idx - 1

        def _parse_cell_ref(cell_ref: str) -> int:
            letters = []
            for ch in cell_ref:
                if "A" <= ch <= "Z":
                    letters.append(ch)
                else:
                    break
            return _col_letters_to_index("".join(letters))

        def _excel_serial_to_date_string(value: float, date1904: bool) -> str:
            if date1904:
                origin = datetime.datetime(1904, 1, 1)
                dt = origin + datetime.timedelta(days=value)
            else:
                origin = datetime.datetime(1899, 12, 30)  # Excel 1900 system offset
                dt = origin + datetime.timedelta(days=value)
            return dt.strftime("%Y-%m-%d")

        # check extension and working directory, then open the file
        file = self._extension_handler(xlsx_file, ".xlsx")
        if self._in_dir and os.path.dirname(file) == "":
            xlsx_file = os.path.join(self._in_dir, file)

        with zipfile.ZipFile(xlsx_file, "r") as z:
            # Load shared strings
            shared = []
            if "xl/sharedStrings.xml" in z.namelist():
                sst = ET.fromstring(z.read("xl/sharedStrings.xml"))
                for si in sst.findall("main:si", NS):
                    t = si.find("main:t", NS)
                    if t is not None and t.text is not None:
                        shared.append(t.text)
                    else:
                        rt = "".join((r.find("main:t", NS).text or "") for r in si.findall("main:r", NS))
                        shared.append(rt)

            # Workbook and sheet resolution
            wb = ET.fromstring(z.read("xl/workbook.xml"))
            date1904 = False
            wbpr = wb.find("main:workbookPr", NS)
            if wbpr is not None and wbpr.attrib.get("date1904") in ("1", "true", "True"):
                date1904 = True

            sheets = wb.find("main:sheets", NS)
            target_rel_id = None
            for s in sheets.findall("main:sheet", NS):
                if s.attrib.get("name") == sheet_name:
                    target_rel_id = s.attrib.get(f"{{{NS['rel']}}}id")
                    break
            if not target_rel_id:
                raise ValueError(f"Sheet '{sheet_name}' not found")

            wb_rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
            target = None
            for r in wb_rels.findall("pkgrel:Relationship", NS):
                if r.attrib.get("Id") == target_rel_id:
                    target = r.attrib["Target"]
                    break
            if not target:
                raise ValueError(f"Could not resolve sheet '{sheet_name}' target")

            sheet_xml = ET.fromstring(z.read(f"xl/{target}"))

            # Styles (for date detection)
            styles_map = {}
            if "xl/styles.xml" in z.namelist():
                styles = ET.fromstring(z.read("xl/styles.xml"))
                numfmt_by_id = {}
                numFmts = styles.find("main:numFmts", NS)
                if numFmts is not None:
                    for nf in numFmts.findall("main:numFmt", NS):
                        nid = nf.attrib.get("numFmtId")
                        code = nf.attrib.get("formatCode", "")
                        if nid:
                            numfmt_by_id[int(nid)] = code
                builtin_date_ids = {14, 15, 16, 17, 22, 27, 30, 36, 45, 46, 47, 57, 58}
                def looks_like_date(code: str) -> bool:
                    c = code.lower()
                    return any(k in c for k in ("y", "m", "d"))
                cellXfs = styles.find("main:cellXfs", NS)
                if cellXfs is not None:
                    for i, xf in enumerate(cellXfs.findall("main:xf", NS)):
                        numFmtId = xf.attrib.get("numFmtId")
                        is_date = False
                        if numFmtId is not None:
                            nid = int(numFmtId)
                            if nid in numfmt_by_id:
                                if looks_like_date(numfmt_by_id[nid]):
                                    is_date = True
                            elif nid in builtin_date_ids:
                                is_date = True
                        styles_map[i] = is_date

            # Parse rows
            rows = []
            sheetData = sheet_xml.find("main:sheetData", NS)
            if sheetData is None:
                rows = []
            else:
                for row in sheetData.findall("main:row", NS):
                    cells = row.findall("main:c", NS)
                    if not cells:
                        rows.append([])
                        continue

                    max_col = -1
                    refs = []
                    for c in cells:
                        cref = c.attrib.get("r")
                        ci = _parse_cell_ref(cref) if cref else len(refs)
                        refs.append((ci, c))
                        if ci > max_col:
                            max_col = ci

                    row_vals = [""] * (max_col + 1)

                    for ci, c in refs:
                        t_attr = c.attrib.get("t")
                        s_attr = c.attrib.get("s")
                        v = c.find("main:v", NS)

                        if t_attr == "inlineStr":
                            is_el = c.find("main:is", NS)
                            t_el = is_el.find("main:t", NS) if is_el is not None else None
                            row_vals[ci] = t_el.text if t_el is not None else ""
                            continue

                        if v is None or v.text is None:
                            row_vals[ci] = ""
                            continue

                        if t_attr == "s":
                            idx = int(v.text)
                            row_vals[ci] = shared[idx] if 0 <= idx < len(shared) else ""
                        else:
                            raw = v.text
                            if s_attr is not None and int(s_attr) in styles_map and styles_map[int(s_attr)]:
                                try:
                                    serial = float(raw)
                                    row_vals[ci] = _excel_serial_to_date_string(serial, date1904)
                                except ValueError:
                                    row_vals[ci] = raw
                            else:
                                row_vals[ci] = raw  # always string otherwise

                    rows.append(row_vals)

        # create table object
        t = self.new_table()

        # headers
        if headers:
            col_headers = rows[0]
            data = [self._row_dtype(row) for row in rows[1:]]
            for h in range(len(col_headers)):  # remove extra quotes from headers
                if col_headers[h].endswith('"') and col_headers[h].startswith('"'):
                    col_headers[h] = col_headers[h][1:-1]
            t._load_on_read(data)
            t._headers = col_headers
        else:
            data = [self._row_dtype(row) for row in rows]
            t._load_on_read(data)
            t.autogenerate_headers()

        t.detect_dtypes()
        return t

    # ============================================================
    # WRITERS
    # ============================================================

    def write_text(self, table, filename, delimiter=","):
        """Attempt to write the data stored in a Table object to a text file ('filename') with a ".txt" extension.
        
        Parameters:
        
        table (object): the Table object to be written.
        
        filename (string): the name of the file to be written.
        
        delimiter (string): the character that separates values in the output file. A comma is the default delimiter.
        """
        data = table.get_data()
        if table.get_headers():
            data.insert(0, table.get_headers())

        if not data:
            print("Table object has no data to be written.")
            return

        # normalize extension
        filename = self._extension_handler(filename, ".txt")

        # Apply output directory
        if self._out_dir and os.path.dirname(filename) == "":
            filename = os.path.join(self._out_dir, filename)

        # Overwrite protection
        if not self._allow_overwrite and self.file_exists(filename):
            raise FileExistsError(f"File '{filename}' already exists and overwrite is disabled. To allow overwriting, set the 'allow_overwrite' attribute to True.")

        # Build all rows as strings in one pass
        lines = [delimiter.join(self._escape_field(val,delimiter) for val in row) + "\n" for row in data]

        # Write all rows at once
        with open(filename, "w", encoding="utf-8") as f:
            f.writelines(lines)

    def write_other_text(self,table,filename,extension,delimiter = ","):
        """Attempt to write the data stored in a Table object to a text file ('filename') with an extension other than ".txt".
        
        Parameters:

        table (object): the Table object to be written.
        
        filename (string): the name of the file to be written.
        
        delimiter (string): the character that separates values in the output file. A comma is the default delimiter.
        
        extension (string): the extension of the file to be written."""
        data = table.get_data()
        if table.get_headers():
            data.insert(0, table.get_headers())

        if not data:
            print("Table object has no data to be written.")
            return

        # normalize extension
        filename = self._extension_handler(filename, target_ext=extension)

        # Apply output directory
        if self._out_dir and os.path.dirname(filename) == "":
            filename = os.path.join(self._out_dir, filename)
        
        # Overwrite protection
        if not self._allow_overwrite and self.file_exists(filename):
            raise FileExistsError(f"File '{filename}' already exists and overwrite is disabled. To allow overwriting, set the 'allow_overwrite' attribute to True.")

        # Build all rows as strings in one pass
        lines = [delimiter.join(self._escape_field(val,delimiter) for val in row) + "\n" for row in data]

        # Write all rows at once
        with open(filename, "w", encoding="utf-8") as f:
            f.writelines(lines)

    def write_csv(self,table,filename):
        """Attempt to write the data stored in a Table object to a CSV file ('filename') with a ".csv" extension.

        Parameters:
        
        table (object): the Table object to be written.
        
        filename (string): the name of the file to be written."""
        data = table.get_data()
        if table.get_headers():
            data.insert(0, table.get_headers())

        if not data:
            print("Table object has no data to be written.")
            return

        # normalize extension
        filename = self._extension_handler(filename, ".csv")

        # Apply output directory
        if self._out_dir and os.path.dirname(filename) == "":
            filename = os.path.join(self._out_dir, filename)
        
        # Overwrite protection
        if not self._allow_overwrite and self.file_exists(filename):
            raise FileExistsError(f"File '{filename}' already exists and overwrite is disabled. To allow overwriting, set the 'allow_overwrite' attribute to True.")

        # Build all rows as strings in one pass
        lines = [','.join(self._escape_field(val,",") for val in row) + "\n" for row in data]

        # Write all rows at once
        with open(filename, "w", encoding="utf-8") as f:
            f.writelines(lines)

    def write_tsv(self, table, filename):
        """Attempt to write the data stored in a Table object to a TSV file ('filename') with a ".tsv" extension.

        Parameters:

        table (object): the Table object to be written.

        filename (string): the name of the file to be written."""
        # Extract data and prepend headers
        data = table.get_data()
        if table.get_headers():
            data.insert(0, table.get_headers())

        if not data:
            print("Table object has no data to be written.")
            return

        # Normalize extension
        filename = self._extension_handler(filename, ".tsv")

        # Apply output directory
        if self._out_dir and os.path.dirname(filename) == "":
            filename = os.path.join(self._out_dir, filename)

        # Overwrite protection
        if not self._allow_overwrite and self.file_exists(filename):
            raise FileExistsError(
                f"File '{filename}' already exists and overwrite is disabled. "
                "To allow overwriting, set the 'allow_overwrite' attribute to True."
            )

        # Build lines
        lines = ["\t".join(self._escape_field(val,"\t") for val in row) + "\n" for row in data]

        # Write file
        with open(filename, "w", encoding="utf-8") as f:
            f.writelines(lines)

    def write_psv(self, table, filename):
        """Attempt to write the data stored in a Table object to a PSV file ('filename') with a ".psv" extension.

        Parameters:

        table (object): the Table object to be written.

        filename (string): the name of the file to be written."""
        # Extract data and prepend headers
        data = table.get_data()
        if table.get_headers():
            data.insert(0, table.get_headers())

        if not data:
            print("Table object has no data to be written.")
            return

        # Normalize extension
        filename = self._extension_handler(filename, ".psv")

        # Apply output directory
        if self._out_dir and os.path.dirname(filename) == "":
            filename = os.path.join(self._out_dir, filename)

        # Overwrite protection
        if not self._allow_overwrite and self.file_exists(filename):
            raise FileExistsError(
                f"File '{filename}' already exists and overwrite is disabled. "
                "To allow overwriting, set the 'allow_overwrite' attribute to True."
            )

        # Build lines
        lines = ["|".join(self._escape_field(val,"|") for val in row) + "\n" for row in data]

        # Write file
        with open(filename, "w", encoding="utf-8") as f:
            f.writelines(lines)

    def write_fixedwidth(self, table, filename, extension = ".txt", widths=None, overflow=None):
        """Attempt to write the data stored in a Table object to a fixed-width text file ('filename') with a ".txt" extension.

        Parameters:

        table (object): the Table object to be written.

        filename (string): the name of the file to be written.

        extension (string): the extension of the file to be read.

        widths (list): a list of character widths of each column. If None, widths will be auto-detected from the data.

        overflow (string): how column overflow should be managed if a value exceeds its specified width. If None, values are silently truncated. If "auto-expand", columns are expanded as needed to fit values. If "strict", an error is raised."""
        # Extract data and prepend headers
        data = table.get_data()
        if table.get_headers():
            data.insert(0, table.get_headers())

        if not data:
            print("Table object has no data to be written.")
            return

        # Normalize extension
        filename = self._extension_handler(filename, extension)

        # Apply output directory
        if self._out_dir and os.path.dirname(filename) == "":
            filename = os.path.join(self._out_dir, filename)

        # Overwrite protection
        if not self._allow_overwrite and self.file_exists(filename):
            raise FileExistsError(
                f"File '{filename}' already exists and overwrite is disabled. "
                "To allow overwriting, set the 'allow_overwrite' attribute to True."
            )

        # AUTO-DETECT WIDTHS IF NOT PROVIDED
        if widths is None:
            # Determine max width of each column
            num_cols = len(data[0])
            widths = [0] * num_cols

            for row in data:
                for i, val in enumerate(row):
                    val_str = str(val)
                    if len(val_str) > widths[i]:
                        widths[i] = len(val_str)

        # HANDLE OVERFLOW BEHAVIOR
        if overflow == "auto-expand":
            # Expand widths to fit the longest value in each column
            num_cols = len(widths)
            for row in data:
                for i in range(num_cols):
                    val_str = str(row[i])
                    if len(val_str) > widths[i]:
                        widths[i] = len(val_str)

        elif overflow == "strict":
            # Raise error if any value exceeds its width
            for r, row in enumerate(data):
                for c, val in enumerate(row):
                    val_str = str(val)
                    if len(val_str) > widths[c]:
                        raise ValueError(
                            f"Fixed-width overflow detected in column {c+1}, row {r+1}: "
                            f"'{val_str}' exceeds width {widths[c]}"
                        )

        # BUILD FIXED-WIDTH ROWS
        lines = []
        for row in data:
            cells = []
            for i, val in enumerate(row):
                val_str = str(val)

                # Truncate if needed (default mode)
                if len(val_str) > widths[i]:
                    val_str = val_str[:widths[i]]

                # Pad right to width
                padded = val_str.ljust(widths[i])
                cells.append(padded)

            # Join cells with no delimiter (fixed-width)
            lines.append("".join(cells) + "\n")

        # WRITE FILE
        with open(filename, "w", encoding="utf-8") as f:
            f.writelines(lines)

    def write_markdown(self, table, filename):
        """Attempt to write the data stored in a Table object to a Markdown file ('filename') with a ".md" extension.

        Parameters:

        table (object): the Table object to be written.

        filename (string): the name of the file to be written."""
        # Extract data and prepend headers
        data = table.get_data()
        if table.get_headers():
            data.insert(0, table.get_headers())

        if not data:
            print("Table object has no data to be written.")
            return

        # Normalize extension
        filename = self._extension_handler(filename, ".md", valid_exts=[".md", ".markdown"])

        # Apply output directory
        if self._out_dir and os.path.dirname(filename) == "":
            filename = os.path.join(self._out_dir, filename)

        # Overwrite protection
        if not self._allow_overwrite and self.file_exists(filename):
            raise FileExistsError(
                f"File '{filename}' already exists and overwrite is disabled. "
                "To allow overwriting, set the 'allow_overwrite' attribute to True."
            )

        # Convert all values to strings
        str_rows = [[str(val) for val in row] for row in data]

        # Determine column count
        num_cols = len(str_rows[0])

        # Header row
        header_row = "| " + " | ".join(str_rows[0]) + " |"

        # Separator row (GitHub-style)
        separator_row = "| " + " | ".join(["---"] * num_cols) + " |"

        # Data rows
        data_rows = [
            "| " + " | ".join(row) + " |"
            for row in str_rows[1:]
        ]

        # Combine all rows
        md_output = "\n".join([header_row, separator_row] + data_rows) + "\n"

        # WRITE FILE
        with open(filename, "w", encoding="utf-8") as f:
            f.write(md_output)

    def write_latex(self, table, filename):
        """Attempt to write the data stored in a Table object to a LaTeX file ('filename') with a ".tex" extension.

        Parameters:

        table (object): the Table object to be written.

        filename (string): the name of the file to be written."""
        # Extract data and prepend headers
        data = table.get_data()
        if table.get_headers():
            data.insert(0, table.get_headers())

        if not data:
            print("Table object has no data to be written.")
            return

        # Normalize extension
        filename = self._extension_handler(filename, ".tex")

        # Apply output directory
        if self._out_dir and os.path.dirname(filename) == "":
            filename = os.path.join(self._out_dir, filename)

        # Overwrite protection
        if not self._allow_overwrite and self.file_exists(filename):
            raise FileExistsError(
                f"File '{filename}' already exists and overwrite is disabled. "
                "To allow overwriting, set the 'allow_overwrite' attribute to True."
            )

        # Convert all values to strings
        str_rows = [[str(val) for val in row] for row in data]

        # Determine number of columns
        num_cols = len(str_rows[0])

        # Column alignment: center all columns (matches common usage)
        col_format = "c" * num_cols

        # Begin environment
        lines = []
        lines.append(r"\begin{tabular}{" + col_format + "}")
        lines.append(r"\hline")

        # Header row
        header_row = " & ".join(str_rows[0]) + r" \\"
        lines.append(header_row)
        lines.append(r"\hline")

        # Data rows
        for row in str_rows[1:]:
            row_tex = " & ".join(row) + r" \\"
            lines.append(row_tex)

        # End environment
        lines.append(r"\hline")
        lines.append(r"\end{tabular}")
        lines.append("")  # final newline

        latex_output = "\n".join(lines)

        with open(filename, "w", encoding="utf-8") as f:
            f.write(latex_output)

    def write_yaml(self, table, filename):
        """Attempt to write the data stored in a Table object to a YAML file ('filename') with a ".yaml" extension. The YAML document will contain a top-level list of dictionaries, where each dictionary represents a row and keys correspond to column headers.

        Parameters:

        table (object): the Table object to be written.

        filename (string): the name of the file to be written."""
        # Extract data and prepend headers
        data = table.get_data()
        headers = table.get_headers()
        if headers:
            data.insert(0, headers)

        if not data:
            print("Table object has no data to be written.")
            return

        # Normalize extension (.yaml preferred, .yml allowed)
        filename = self._extension_handler(filename, ".yaml", valid_exts=[".yaml", ".yml"])

        # Apply output directory
        if self._out_dir and os.path.dirname(filename) == "":
            filename = os.path.join(self._out_dir, filename)

        # Overwrite protection
        if not self._allow_overwrite and self.file_exists(filename):
            raise FileExistsError(
                f"File '{filename}' already exists and overwrite is disabled. "
                "To allow overwriting, set the 'allow_overwrite' attribute to True."
            )

        # First row is headers
        col_headers = data[0]
        body_rows = data[1:]

        yaml_lines = []

        for row in body_rows:
            yaml_lines.append("-")  # start of dict
            for key, val in zip(col_headers, row):
                yaml_lines.append(f"  {key}: {val}")

        yaml_output = "\n".join(yaml_lines) + "\n"

        # WRITE FILE
        with open(filename, "w", encoding="utf-8") as f:
            f.write(yaml_output)

    def write_ini(self, table, filename):
        """Attempt to write the data stored in a Table object to an INI-style configuration file ('filename') with a ".ini" extension. Each row in the Table becomes a section, and each column becomes a key-value pair within that section.

        Parameters:

        table (object): the Table object to be written.

        filename (string): the name of the file to be written."""
        
        # Extract data and prepend headers
        data = table.get_data()
        headers = table.get_headers()
        if headers:
            data.insert(0, headers)

        if not data:
            print("Table object has no data to be written.")
            return

        # Normalize extension
        filename = self._extension_handler(filename, ".ini", valid_exts=[".ini", ".cfg", ".conf"])

        # Apply output directory
        if self._out_dir and os.path.dirname(filename) == "":
            filename = os.path.join(self._out_dir, filename)

        # Overwrite protection
        if not self._allow_overwrite and self.file_exists(filename):
            raise FileExistsError(
                f"File '{filename}' already exists and overwrite is disabled. "
                "To allow overwriting, set the 'allow_overwrite' attribute to True."
            )

        # First row is headers
        col_headers = data[0]
        body_rows = data[1:]

        ini_lines = []

        for idx, row in enumerate(body_rows):
            # Section name: Section1, Section2, ...
            section_name = f"Section{idx+1}"
            ini_lines.append(f"[{section_name}]")

            # Key-value pairs
            for key, val in zip(col_headers, row):
                ini_lines.append(f"{key} = {val}")

            ini_lines.append("")  # blank line between sections

        ini_output = "\n".join(ini_lines)

        # WRITE FILE
        with open(filename, "w", encoding="utf-8") as f:
            f.write(ini_output)

    def write_xml(self, table, filename):
        """Attempt to write the data stored in a Table object to an XML file ('filename') with a ".xml" extension. The XML document will contain a single <table> element, with each row represented as a <row> element containing child tags corresponding to column headers.

        Parameters:

        table (object): the Table object to be written.

        filename (string): the name of the file to be written."""
        # Extract data and prepend headers
        data = table.get_data()
        headers = table.get_headers()
        if headers:
            data.insert(0, headers)

        if not data:
            print("Table object has no data to be written.")
            return

        # Normalize extension
        filename = self._extension_handler(filename, ".xml")

        # Apply output directory
        if self._out_dir and os.path.dirname(filename) == "":
            filename = os.path.join(self._out_dir, filename)

        # Overwrite protection
        if not self._allow_overwrite and self.file_exists(filename):
            raise FileExistsError(
                f"File '{filename}' already exists and overwrite is disabled. "
                "To allow overwriting, set the 'allow_overwrite' attribute to True."
            )

        # First row is headers
        col_headers = data[0]
        body_rows = data[1:]

        xml_lines = []
        xml_lines.append("<table>")

        for row in body_rows:
            xml_lines.append("  <row>")
            for key, val in zip(col_headers, row):
                xml_lines.append(f"    <{key}>{val}</{key}>")
            xml_lines.append("  </row>")

        xml_lines.append("</table>")
        xml_lines.append("")  # final newline

        xml_output = "\n".join(xml_lines)

        # WRITE FILE
        with open(filename, "w", encoding="utf-8") as f:
            f.write(xml_output)

    def write_json(self,table,filename):
        """Attempt to write the data stored in a Table object to a JSON file ('filename') with a ".json" extension. This JSON file will have a single object. The object's keys will be the Table column headers, and the object's values will be arrays of column values.
        
        Parameters:
        
        table (object): the Table object to be written.
        
        filename (string): the name of the file to be written."""
        import json
        if table.is_empty():
            print("Table object has no data to be written.")
            return

        # normalize extension
        filename = self._extension_handler(filename, ".json")

        # Apply output directory
        if self._out_dir and os.path.dirname(filename) == "":
            filename = os.path.join(self._out_dir, filename)
        
        # Overwrite protection
        if not self._allow_overwrite and self.file_exists(filename):
            raise FileExistsError(f"File '{filename}' already exists and overwrite is disabled. To allow overwriting, set the 'allow_overwrite' attribute to True.")

        # convert table to dict
        jdict = self.table_to_dict(table)

        # safe file writing
        with open(filename, "w", encoding="utf-8") as f:
            json.dump(jdict, f, ensure_ascii=False, indent=4)

    def write_html(self,table,filename,style = "standard"):
        """Attempt to write the data stored in a Table object to a table in an HTML file ('filename') with the ".html" extension.
        
        Parameters:

        table (object): the Table object to be written.
        
        filename (string): the name of the file to be written.
        
        style (string): the style of the output table. The style may be specified as "standard" to enclose all cells in horizontal and vertical borders, or "scientific" to remove vertical borders from cells."""
        data = table.get_data()
        if table.get_headers():
            data.insert(0, table.get_headers())

        if not data:
            print("Table object has no data to be written.")
            return

        # normalize extension
        filename = self._extension_handler(filename, ".html", valid_exts=[".html", ".htm"])

        # Apply output directory
        if self._out_dir and os.path.dirname(filename) == "":
            filename = os.path.join(self._out_dir, filename)
        
        # Overwrite protection
        if not self._allow_overwrite and self.file_exists(filename):
            raise FileExistsError(f"File '{filename}' already exists and overwrite is disabled. To allow overwriting, set the 'allow_overwrite' attribute to True.")

        # choose style
        if style == "scientific":
            header = """<!doctype html>
    <html lang="en">
    <head>
        <meta charset="utf-8">
        <title>TableTools Table</title>
        <style type="text/css">
            table {border-collapse: collapse; border-left: none; border-right: none; font-size: 20px;}
            td, th {border: 1px solid #000000; border-left: none; border-right: none; text-align: center; padding: 8px;}
        </style>
    </head>
    <body>
        <table align="left">"""
        else:  # default standard
            header = """<!doctype html>
    <html lang="en">
    <head>
        <meta charset="utf-8">
        <title>TableTools Table</title>
        <style type="text/css">
            table {border-collapse: collapse; font-size: 20px;}
            td, th {border: 1px solid #000000; text-align: center; padding: 4px;}
        </style>
    </head>
    <body>
        <table align="left">"""

        rows_html = []
        for i, row in enumerate(data):
            tag = "th" if i == 0 else "td"
            cells = "".join(f"<{tag}>{str(val)}</{tag}>" for val in row)
            rows_html.append(f"<tr>{cells}</tr>")

        footer = "</table>\n</body>\n</html>"
        html_output = "\n".join([header] + rows_html + [footer])

        with open(filename, "w", encoding="utf-8") as f:
            f.write(html_output)

    def write_r_dput(self, table, filename):
        """Attempt to write the data stored in a Table object to an R dput() representation of a data frame ('filename') with a ".r" extension. The output will be a structure(list(...)) object where each column is represented as a named vector using R's c(...) syntax.

        Parameters:

        table (object): the Table object to be written.

        filename (string): the name of the file to be written."""
        # Extract data and prepend headers
        data = table.get_data()
        headers = table.get_headers()
        if headers:
            data.insert(0, headers)

        if not data:
            print("Table object has no data to be written.")
            return

        # Normalize extension
        filename = self._extension_handler(filename, ".r", valid_exts=[".r", ".dput"])

        # Apply output directory
        if self._out_dir and os.path.dirname(filename) == "":
            filename = os.path.join(self._out_dir, filename)

        # Overwrite protection
        if not self._allow_overwrite and self.file_exists(filename):
            raise FileExistsError(
                f"File '{filename}' already exists and overwrite is disabled. "
                "To allow overwriting, set the 'allow_overwrite' attribute to True."
            )

        # First row is headers
        col_headers = data[0]
        body_rows = data[1:]

        # Build column-wise vectors
        columns = {h: [] for h in col_headers}

        for row in body_rows:
            for h, val in zip(col_headers, row):
                columns[h].append(val)

        # Sort column names (reader sorts them)
        sorted_headers = sorted(col_headers)

        # Format values according to R rules
        def r_format(value):
            if value is None:
                return "NA"
            if isinstance(value, bool):
                return "TRUE" if value else "FALSE"
            if isinstance(value, (int, float)):
                return str(value)
            # Strings must be quoted
            return f"\"{value}\""

        # Build list(...) entries
        list_entries = []
        for h in sorted_headers:
            vec = columns[h]
            formatted = ", ".join(r_format(v) for v in vec)
            entry = f"{h} = c({formatted})"
            list_entries.append(entry)

        # Join entries with commas
        list_body = ",\n    ".join(list_entries)

        # Wrap in structure(list(...))
        r_output = "structure(list(\n    " + list_body + "\n))\n"

        # WRITE FILE
        with open(filename, "w", encoding="utf-8") as f:
            f.write(r_output)

    def write_sql(self,table,filename,table_name):
        """Attempt to write the data stored in a Table object to a table in a new SQL database file ('filename'), or add a new table to an existing SQL database file. Data headers should not have spaces. Spaces in the data headers may be corrected with the remove_header_space function in the Table object.
        
        Parameters:

        table (object): the Table object to be written.

        filename (string): the name of the database file to be written. If the specified filename points to an existing SQL database file, the function will attempt to write the data to the existing database file as a new table. Possible database file extensions include ".db", ".db2", ".db3", ".sqlite", ".sqlite2", and ".sqlite3", but if a database extension is not included in the file name, the file will be written as a ".sqlite" file.
        
        table_name (string): the name of the table to be written to the database file. The table name may have underscores but should not have spaces. If the table name already exists in the database, the table will be overwritten, but other tables in the database will remain intact."""
        data = table.get_data()
        if table.get_headers():
            data.insert(0, table.get_headers())
        dtypes = table._dtypes[:]
        if not data:
            print("Table object has no data to be written.")
        else:
            import sqlite3
            # normalize extension
            filename = self._extension_handler( filename, target_ext=".sqlite", valid_exts=[".db", ".db2", ".db3", ".sqlite", ".sqlite2", ".sqlite3"])

            # Apply output directory
            if self._out_dir and os.path.dirname(filename) == "":
                filename = os.path.join(self._out_dir, filename)

            headers = data[0]

            # map dtypes to SQLite types
            dtype_map = {"integer": "INTEGER", "float": "REAL", "string": "TEXT"}
            h_types = [f'"{headers[i]}" {dtype_map.get(dtypes[i], "TEXT")}' for i in range(len(headers))]

            qs = ",".join(["?"] * len(headers))

            # Overwrite protection
            if not self._allow_overwrite and self.file_exists(filename):
                raise FileExistsError(f"File '{filename}' already exists and overwrite is disabled. To allow overwriting, set the 'allow_overwrite' attribute to True.")

            # safe connection handling
            with sqlite3.connect(filename) as conn:
                crs = conn.cursor()
                try:
                    crs.execute(f"CREATE TABLE {table_name} ({','.join(h_types)})")
                    crs.executemany(f"INSERT INTO {table_name} VALUES ({qs})", data[1:])
                except sqlite3.OperationalError:
                    # overwrite existing table
                    crs.execute(f"DROP TABLE {table_name}")
                    crs.execute(f"CREATE TABLE {table_name} ({','.join(h_types)})")
                    crs.executemany(f"INSERT INTO {table_name} VALUES ({qs})", data[1:])
                conn.commit()
   
    def write_dbf(self,table,filename, csize, cdeci):
        """Attempt to write the data stored in a Table object to a Dbase/Xbase file ('filename') with a ".dbf" extension. Headers should not be longer than 10 characters. Headers longer than 10 characters will be reduced to 10 characters. Empty strings '' will be considered Null values.
        
        Parameters:

        table (object): the Table object to be written.
        
        filename (string): the name of the file to be written. Leave this parameter blank to overwrite the original file, if specified.
        
        csize (list): a list of the sizes (the maximum length of values) of each column, in the order the columns appear in the data.
        
        cdeci (list): a list of the number of decimal places of each column, in the order the columns appear in the data. Columns of float values should have a value greater than 0 to indicate the number of decimals. All other columns should be given a value of 0."""
        # The core component of this Dbase/Xbase writer (packing the data) is based on the works of:
        # Raymond Hettinger (https://code.activestate.com/recipes/362715-dbf/)
        # Tomas Nordin (https://code.activestate.com/recipes/580696-dbf-reader-and-writer-selective-fields-and-nullrep/)

        data = table.get_data()
        if table.get_headers():
            data.insert(0, table.get_headers())
        dtypes = table._dtypes[:]
        if not data:
            print("Table object has no data to be written.")
        else:
            import struct
            
            # normalize extension
            filename = self._extension_handler(filename, ".dbf")

            # Apply output directory
            if self._out_dir and os.path.dirname(filename) == "":
                filename = os.path.join(self._out_dir, filename)

            # map dtypes to DBF field types
            dtype_map = {"integer": "N", "float": "N", "string": "C"}
            dtypes = [dtype_map.get(dt, "C") for dt in dtypes]
            field_specs = list(zip(dtypes, csize, cdeci))

            headers = data[0]
            num_rows = len(data) - 1
            num_cols = len(headers)
            lenheader = num_cols * 32 + 33
            len_record = sum(size for _, size, _ in field_specs) + 1

            # file header
            now = datetime.datetime.now()
            ver = 3
            yr, mon, day = now.year - 1900, now.month, now.day
            file_header = struct.pack("<BBBBLHH20x", ver, yr, mon, day, num_rows, lenheader, len_record)

            # Overwrite protection
            if not self._allow_overwrite and self.file_exists(filename):
                raise FileExistsError(f"File '{filename}' already exists and overwrite is disabled. To allow overwriting, set the 'allow_overwrite' attribute to True.")

            with open(filename, "wb") as f:
                f.write(file_header)

                # field descriptors
                for name, (typ, size, deci) in zip(headers, field_specs):
                    name = name[:10].ljust(11, "\x00")  # truncate to 10 chars
                    f.write(struct.pack("<11sc4xBB14x", name.encode("ascii", errors="ignore"),
                                        typ.encode("ascii"), size, deci))

                # header terminator
                f.write(b"\r")

                # records
                for row in data[1:]:
                    f.write(b" ")  # deletion flag
                    for (typ, size, deci), value in zip(field_specs, row):
                        if typ == "N":
                            val = str(value).rjust(size, " ")
                        elif typ == "D" and hasattr(value, "strftime"):
                            val = value.strftime("%Y%m%d")
                        elif typ == "L":
                            val = str(value)[0].upper()
                        else:
                            val = str(value)[:size].ljust(size, " ")
                        f.write(val.encode("ascii", errors="ignore"))

                # EOF marker
                f.write(b"\x1A")
