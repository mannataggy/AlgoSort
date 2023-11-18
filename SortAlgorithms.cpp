#include "SortAlgorithms.h"

#include <atomic>
#include <limits>
#include <tuple>

using SortableIterator = std::vector<Sortable>::iterator;
using std::tuple;

//
// ──────────────────────────────────────────────────────────────────────────────
//   :::::: S O R T   A L G O R I T H M S : :  :   :    :     :        :          
// ───────────────────────────────────────────────────────────────────────────────
//

/**
 *  Bubble sort
 */ 
int algo::bubbleSort(std::vector<Sortable>& sortElements, int timeSleep, const std::atomic<bool>& interrupt) {
	int numOfComparisons = 0;

	for (int n = 0; n < sortElements.size() - 1; n++) {
		if (interrupt) {
			return numOfComparisons;
		}

		if (sortElements[n].value > sortElements[n + 1].value) {
			algoUtils::swap(sortElements, timeSleep, sortElements[n], sortElements[n+1]);
		}
		numOfComparisons++;
	}

	return numOfComparisons;
}

/**
 *  Selection sort
 */
int algo::selectionSort(std::vector<Sortable>& sortElements, int timeSleep, const std::atomic<bool>& interrupt) {
	int numOfComparisons = 0;

	for (int n = 0; n <= sortElements.size() - 1; n++) {
		for (int j = n; j <= sortElements.size() - 1; j++) {
			if (interrupt) {
				return numOfComparisons;
			}

			if (sortElements[n].value > sortElements[j].value) {
				algoUtils::swap(sortElements, timeSleep, sortElements[n], sortElements[j]);

			}
			numOfComparisons++;
		}
	}

	return numOfComparisons;
}

/**
 * Insertion sort
 */
int algo::insertionSort(std::vector<Sortable>& sortElements, int timeSleep, const std::atomic<bool>& interrupt) {
	int numOfComparisons = 0;

	for (int n = 1; n < sortElements.size(); n++)
	{
		auto temp = sortElements[n];
		int j = n - 1;
		while (j >= 0 && temp.value <= sortElements[j].value)
		{
			if (interrupt) {
				return numOfComparisons;
			}

			sortElements[j].color = sf::Color::Red;

			sortElements[j + 1] = sortElements[j];

			sf::sleep(sf::milliseconds(timeSleep));
			sortElements[j+1].color = sf::Color::White;

			j = j - 1;

			numOfComparisons++;

		}
		sortElements[j + 1] = temp;
		numOfComparisons++;
	}

	return numOfComparisons;
}

/**
 * @brief QuickSort Partition step. Iterators follow STL style for ranges.
 *
 * @param parent  The parent array, needed for the swap utilities.
 * @param beg First element of the sub array
 * @param end "One past the last" element of the sub array
 * @param timeSleep pauses the thread for this many ms
 * @param interrupt Bool to stop the sort process
 * @return tuple, first is the pivot, second is the number of comparisons
 */
static int quickSortHelper(std::vector<Sortable>& parent, SortableIterator beg, SortableIterator end, int timeSleep, const std::atomic<bool>& interrupt);
static tuple<SortableIterator, int> quickSortPartition(std::vector<Sortable>& parent, SortableIterator beg, SortableIterator end, int timeSleep, const std::atomic<bool>& interrupt);

static int quickSortHelper(std::vector<Sortable>& parent, SortableIterator beg, SortableIterator end, int timeSleep, const std::atomic<bool>& interrupt) {
	if (interrupt) {
		return 0;
	}

	// base case
	if (end - beg < 2) return 0;

	SortableIterator pivot;
	int numOfComparisons;
	std::tie(pivot, numOfComparisons) = quickSortPartition(parent, beg, end, timeSleep, interrupt);

	return numOfComparisons +
		quickSortHelper(parent, beg, pivot, timeSleep, interrupt) +
		quickSortHelper(parent, pivot + 1, end, timeSleep, interrupt);

}

static tuple<SortableIterator, int> quickSortPartition( std::vector<Sortable>& parent, SortableIterator beg, SortableIterator end, int timeSleep, const std::atomic<bool>& interrupt) {
	auto pivot = end - 1;
	int numOfComparisons = 0;

	auto lhs = beg;
	for (auto rhs = lhs; rhs != pivot; ++rhs) {
		if (interrupt) {
			return std::make_tuple(lhs, numOfComparisons);
		}

		++numOfComparisons;
		if (rhs->value <= pivot->value) {
			algoUtils::swap(parent, timeSleep, *lhs, *rhs);
			++lhs;
		}
	}
	algoUtils::swap(parent, timeSleep, *lhs, *pivot);
	return tuple<SortableIterator, int>{lhs, numOfComparisons};
}

int algo::quickSort(std::vector<Sortable>& sortElements, int timeSleep, const std::atomic<bool>& interrupt) {
	return quickSortHelper(sortElements, sortElements.begin(), sortElements.end(), timeSleep, interrupt);
}

/**
 * @brief Heapify the heap rooted at index 'i'.
 * @param n Size of the heap
 * @param i Root index of the heap
 */
static void heapify(std::vector<Sortable>& sortElements, int timeSleep, const std::atomic<bool>& interrupt, int n, int i, int& numOfComparisons) {
    int largest = i; // Initialize largest as root
    int left = 2 * i + 1; // left = 2*i + 1
    int right = 2 * i + 2; // right = 2*i + 2

    // If left child is larger than root
    if (left < n && sortElements[left].value > sortElements[largest].value) {
        largest = left;
    }

    // If right child is larger than largest so far
    if (right < n && sortElements[right].value > sortElements[largest].value) {
        largest = right;
    }

    // If largest is not root
    if (largest != i) {
        algoUtils::swap(sortElements, timeSleep, sortElements[i], sortElements[largest]);
        numOfComparisons++;

        // Recursively heapify the affected sub-tree
        heapify(sortElements, timeSleep, interrupt, n, largest, numOfComparisons);
    }
}

/**
 * @brief Heap Sort algorithm
 */
int algo::heapSort(std::vector<Sortable>& sortElements, int timeSleep, const std::atomic<bool>& interrupt) {
    int numOfComparisons = 0;

    // Build heap (rearrange array)
    int n = sortElements.size();
    for (int i = n / 2 - 1; i >= 0; i--) {
        heapify(sortElements, timeSleep, interrupt, n, i, numOfComparisons);
    }

    // One by one extract an element from heap
    for (int i = n - 1; i > 0; i--) {
        // Move current root to end
        algoUtils::swap(sortElements, timeSleep, sortElements[0], sortElements[i]);
        numOfComparisons++;

        // call max heapify on the reduced heap
        heapify(sortElements, timeSleep, interrupt, i, 0, numOfComparisons);
    }

    return numOfComparisons;
}

/**
 * @brief Merge two subarrays of sortElements[].
 * @param l Left index of the subarray
 * @param m Middle index of the subarray
 * @param r Right index of the subarray
 */
static void merge(std::vector<Sortable>& sortElements, int timeSleep, const std::atomic<bool>& interrupt, int l, int m, int r, int& numOfComparisons) {
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    // Create temporary arrays
    std::vector<Sortable> L(n1), R(n2);

    // Copy data to temporary arrays L[] and R[]
    for (i = 0; i < n1; i++) {
        L[i] = sortElements[l + i];
    }
    for (j = 0; j < n2; j++) {
        R[j] = sortElements[m + 1 + j];
    }

    // Merge the temporary arrays back into sortElements[l..r]
    i = 0;
    j = 0;
    k = l;
    while (i < n1 && j < n2 && !interrupt) {
        if (L[i].value <= R[j].value) {
            algoUtils::swap(sortElements, timeSleep, sortElements[k], L[i]);
            i++;
        } else {
            algoUtils::swap(sortElements, timeSleep, sortElements[k], R[j]);
            j++;
        }
        numOfComparisons++;
        k++;
    }

    // Copy the remaining elements of L[], if there are any
    while (i < n1 && !interrupt) {
        algoUtils::swap(sortElements, timeSleep, sortElements[k], L[i]);
        i++;
        k++;
        numOfComparisons++;
    }

    // Copy the remaining elements of R[], if there are any
    while (j < n2 && !interrupt) {
        algoUtils::swap(sortElements, timeSleep, sortElements[k], R[j]);
        j++;
        k++;
        numOfComparisons++;
    }
}

/**
 * @brief Merge Sort algorithm (recursive function)
 * @param l Left index of the subarray
 * @param r Right index of the subarray
 */
static void mergeSortHelper(std::vector<Sortable>& sortElements, int timeSleep, const std::atomic<bool>& interrupt, int l, int r, int& numOfComparisons) {
    if (l < r && !interrupt) {
        // Same as (l+r)/2, but avoids overflow
        int m = l + (r - l) / 2;

        // Sort first and second halves
        mergeSortHelper(sortElements, timeSleep, interrupt, l, m, numOfComparisons);
        mergeSortHelper(sortElements, timeSleep, interrupt, m + 1, r, numOfComparisons);

        // Merge the sorted halves
        merge(sortElements, timeSleep, interrupt, l, m, r, numOfComparisons);
    }
}

/**
 * Merge Sort algorithm
 */
int algo::mergeSort(std::vector<Sortable>& sortElements, int timeSleep, const std::atomic<bool>& interrupt) {
    int numOfComparisons = 0;
    int n = sortElements.size();

    // Start the recursive merge sort
    mergeSortHelper(sortElements, timeSleep, interrupt, 0, n - 1, numOfComparisons);

    return numOfComparisons;
}

//
// ────────────────────────────────────────────────────────────────────
//   :::::: U T I L I T I E S : :  :   :    :     :        :          :
// ────────────────────────────────────────────────────────────────────
//

/**
 * @brief algoUtils::swap - Utility function for swapping elements and changing the colors on swap.
 *
 * @param sortElements The array where the elements being swapped are located
 * @param timeSleep Time in milliseconds to wait between swaps
 * @param el1 First element to be swapped
 * @param el2 Second element to be swapped
 */
void algoUtils::swap(std::vector<Sortable>& sortElements, int timeSleep, Sortable& el1, Sortable& el2) {
	el1.color = sf::Color::Red;
	el2.color = sf::Color::Red;

	auto currElement = el1;
	auto tempElement = el2;
	el1 = tempElement;
	el2 = currElement;

	sf::sleep(sf::milliseconds(timeSleep));

	el1.color = sf::Color::White;
	el2.color = sf::Color::White;
}
