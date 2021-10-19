"""
A collection of sorting algorithms

Included sorting algorithims:
    1. Insertion sort
        - https://en.wikipedia.org/wiki/Insertion_sort
    2. Merge sort
        - https://en.wikipedia.org/wiki/Merge_sort
    3. Bubble sort
        - https://en.wikipedia.org/wiki/Bubble_sort
    4. Max-heap sort
        - https://www.geeksforgeeks.org/python-program-for-heap-sort/
    5. Quick-sort
        - https://www.geeksforgeeks.org/python-program-for-quicksort/

Author: Matthew R. DeVerna
"""
import warnings

from .utils import check_array

def insertion_sort(given_array):
    """
    Implement the **insertion sort** algorithm to sort `given_array` in place.

    Iteratively sort `given_array` moving from left to right. This algorithm
    is an efficient algorithm for sorting a small number of elements.

    Complexity:
    ----------
    - O(n^2)
    Note: Merge sort run time complexity ( O(n*log(n)) ) is beneficial for
        larger arrays. Due to run time constants, insertion sort is likely
        beneficial for shorter arrays.

    Parameters:
    ----------
    - given_array (list) : a numerical sequence

    Exceptions:
    ----------
    - TypeError
    """
    # Ensure array is a list and contains only numeric values
    check_array(given_array)

    for j in range(1, len(given_array)):
        # Variable key represents the element under consideration.
        key = given_array[j]
        i = j - 1
        # While i is not exceeding the left boundary of a list
        # And the element, which has the index of i, is greater than key.
        while i >= 0 and given_array[i] > key:
            # "Shift" the elements upwards
            # as long as the condition remains True.
            # Note, that we are just copying an element with index i
            # to the right by one position.
            given_array[i + 1] = given_array[i]
            i = i - 1
        # Insert key into the proper position.
        given_array[i + 1] = key

def merge_sort(given_array):
    """
    Implement the **merge sort** algorithm to sort `given_array` in place.

    Sort `given_array` **inplace** in ascending order by recursively dividing
    the it into two smaller lists, sorting those, and then merging them all back
    together later.

    Complexity:
    ----------
    - O(n*log(n))
        Note: Although merge sort run time complexity is less than insertion sort
        ( O(n^2) ), the constant factors in insertion sort can make it faster
        in practice for small problem sizes on many machines.

    Parameters:
    ----------
    - given_array (list) : a numerical sequence

    Exceptions:
    ----------
    - TypeError
    """
    # Ensure array is a list and contains only numeric values
    check_array(given_array)

    size = len(given_array)
    if size > 1:

        # Split the array in two
        middle = size // 2
        left_arr = given_array[:middle]
        right_arr = given_array[middle:]

        # Recursively call this function
        merge_sort(left_arr)
        merge_sort(right_arr)

        merged_idx = 0

        # While there are elements in either the left or right
        # arrays we check which is larger and append the smallest
        while (len(left_arr) > 0) or (len(right_arr) > 0):

            # If either list is empty, take from the other
            if len(left_arr) == 0:
                given_array[merged_idx] = right_arr.pop(0)
                merged_idx += 1
            elif len(right_arr) == 0:
                given_array[merged_idx] = left_arr.pop(0)
                merged_idx += 1

            # If both lists have elements...

            # Take right arr[0], if left_arr[0] is larger
            elif left_arr[0] > right_arr[0]:
                given_array[merged_idx] = right_arr.pop(0)
                merged_idx += 1

            # Otherwise, take the first element of left_arr
            else:
                given_array[merged_idx] = left_arr.pop(0)
                merged_idx += 1

def bubble_sort(given_array):
    """
    Implement the **bubble sort** algorithm to sort `given_array` in place in 
    ascending order.

    Repeatedly steps through a list, comparing adjacent elements, and
    swapping them if they are in the wrong order.

    Complexity:
    ----------
    - O(n^2)
        Note: Other О(n^2) sorting algorithms, such as insertion sort,
        generally run faster than bubble sort, and are no more complex.
        **Therefore, bubble sort is not a practical sorting algorithm.**

    Parameters:
    ----------
    - given_array (list) : a numerical sequence

    Exceptions:
    ----------
    - TypeError
    """
    # Ensure array is a list and contains only numeric values
    check_array(given_array)

    n = len(given_array)

    # Check adjacent items and rearrange them if the item on the left is
    # bigger than the item on the right
    for i in range(n-1):
        for j in range(n-i-1):
            if(given_array[j] > given_array[j+1]):
                given_array[j], given_array[j+1] = given_array[j+1], given_array[j]

def _max_heapify(given_array, n, i):
    """
    Transform **inplace** `given_array` into a max-heap array ordering.
        Note: This is only meant to be utilized within `max_heap_sort`.

    Parameters:
    ----------
    - given_array (list) : the array to transform
    - n (int) : the size of the array
    - i (int) : index of largest element (initialized as root (i.e., 0))

    Exceptions:
    ----------
    - TypeError
    """
    # Ensure array is a list and contains only numeric values
    check_array(given_array)

    if not isinstance(n,int):
        raise TypeError(f"`n` must be an integer. Currently it is <{type(n)}>")
    if not isinstance(i,int):
        raise TypeError(f"`i` must be an integer. Currently it is <{type(i)}>")

    largest = i           # Initialize largest as root
    l_idx = 2 * i + 1     # left = 2*i + 1
    r_idx = 2 * i + 2     # right = 2*i + 2

    # If the left child is not the last node and it
    # is > the current largest element --> set new `largest` index
    if (l_idx < n) and (given_array[i] < given_array[l_idx]):
        largest = l_idx

    # If the right child is not at the last node and it
    # is > the current largest element --> set new `largest` index
    if (r_idx < n) and (given_array[largest] < given_array[r_idx]):
        largest = r_idx

    # If we have a new largest node, switch it with the root
    if largest != i:
        given_array[i], given_array[largest] = given_array[largest], given_array[i]

        # Heapify the root.
        _max_heapify(given_array, n, largest)

def max_heap_sort(given_array):
    """
    Sort `given_array` **inplace** in ascending order by first converting it
    to a max-heap structure.

    Complexity:
    ----------
    - O(nlogn)

    Parameters:
    ----------
    - given_array (list) : the array to sort

    Exceptions:
    ----------
    - TypeError
    """
    # Ensure array is a list and contains only numeric values
    check_array(given_array)

    n = len(given_array)

    # Build a maxheap.
    # Since last parent will be at ((n//2)-1) we can start at that location.
    last_parent_idx = (n // 2) - 1
    for i in range(last_parent_idx, -1, -1):
        _max_heapify(given_array, n, i)

    # One by one extract elements
    for i in range(n-1, 0, -1):
        given_array[i], given_array[0] = given_array[0], given_array[i]
        _max_heapify(given_array, i, 0)

def _partition(given_array, low, high):
    """
    A method for partitioning an array (inplace) into two groups.
    The last value of `given_array` is used to put all elements
    of lower (or equal) value on the left of that element and all elements
    of higher value on the right of that element.

    Note: Meant to be utilized within the `quick_sort` function but can
        be used on it's own

    Parameters:
    ----------
    - given_array (list) : the array to be paritioned
    - low (int) : the lowest value index (typically zero)
    - high (int) : the highest value index (typically, len(given_array)-1)

    Exceptions:
    ----------
    - TypeError
    """
    # Ensure array is a list and contains only numeric values
    check_array(given_array)

    if not isinstance(low,int):
        raise TypeError(f"`low` must be an integer. Currently it is <{type(low)}>")
    if not isinstance(high,int):
        raise TypeError(f"`high` must be an integer. Currently it is <{type(high)}>")

    # 1 minus the lowest index to
    #   know where we are in the given_arrayay
    i = low-1

    # Pivot element to compare other elements too
    pivot = given_array[high]

    for j in range(low,high):

        # If jth item is less than or equal to pivot element, swap them
        #   and increment `i` by one
        if given_array[j] <= pivot:
            i += 1
            given_array[j], given_array[i] = given_array[i], given_array[j]

    # Swap the pivot into it's proper position
    #   now that everything is sorted
    given_array[i+1], given_array[high] = given_array[high], given_array[i+1]
    return (i+1)

def quick_sort(given_array, low, high):
    """
    Sort `given_array` **inplace** in ascending order using the quick sort algorithm.

    Complexity:
    ----------
    - Worst: O(n^2)
        - Occurs when the partition process always picks the
        greatest or smallest element as pivot
    - Average: O(nlogn)

    Parameters:
    ----------
    - given_array (list) : the array to be paritioned
    - low (int) : the lowest value index (typically zero)
    - high (int) : the highest value index (typically, len(given_array)-1)

    Exceptions:
    ----------
    - TypeError
    """
    # Ensure array is a list and contains only numeric values
    check_array(given_array)

    if not isinstance(low,int):
        raise TypeError(f"`low` must be an integer. Currently it is <{type(low)}>")
    if not isinstance(high,int):
        raise TypeError(f"`high` must be an integer. Currently it is <{type(high)}>")

    if len(given_array) == 1:
        return given_array
    if low < high:

        # piv is partitioning index, given_array[p] is now
        # at right place
        piv = _partition(given_array, low, high)

        # Separately sort elements before
        # partition and after partition
        quick_sort(given_array, low, piv-1)
        quick_sort(given_array, piv+1, high)

def counting_sort(given_array):
    """
    Sort `given_array` in ascending order via the counting sort algorithm.
    Only accepts positive numbers

    Complexity:
    ----------
    - O(n) - only when the maximum value in `given_array` is < len(given_array)*2

    Parameters:
    ----------
    - given_array : a sequence to be sorted
    """
    # Ensure array is a list and contains only numeric values
    check_array(given_array)

    if any(val<0 for val in given_array):
        raise TypeError("This function does not support negative values.")

    # Set needed parameters and lists
    max_val = max(given_array)
    count_arr = [0] * (max_val+1)
    new_arr = [0] * len(given_array)

    # Throw complexity warning
    if max_val >= len(given_array)*2:
        warnings.warn(
            "\nCOMPLEXITY WARNING:\n"
            "The maximum value within the provided array "
            "is > len(array)*2 so running time will "
            "no longer be O(n)."
        )

    # Here we increment one for given_array[idx] in count_arr[idx]
    for idx in range(0, len(given_array)):
        count_arr[given_array[idx]] = count_arr[given_array[idx]] + 1

    # With the above list, we  determine how many
    # elements are less than or equal idx by keeping a running sum
    for idx in range(1,max_val+1):
        count_arr[idx] = count_arr[idx] + count_arr[idx-1]

    # Place each element given_array[idx] into its correct
    # sorted position in the `new_arr`
    for idx in range(len(given_array)-1,-1,-1):
        new_arr[count_arr[given_array[idx]]-1] = given_array[idx]
        count_arr[given_array[idx]] = count_arr[given_array[idx]]-1

    return new_arr