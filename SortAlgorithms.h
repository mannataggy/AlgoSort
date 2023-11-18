#pragma once
#include <atomic>
#include <vector>
#include "Sortable.h"

namespace algo {
	int bubbleSort(std::vector<Sortable>& sortElements, int timeSleep, const std::atomic<bool>& interrupt);
	int selectionSort(std::vector<Sortable>& sortElements, int timeSleep, const std::atomic<bool>& interrupt);
	int insertionSort(std::vector<Sortable>& sortElements, int timeSleep, const std::atomic<bool>& interrupt);
	int quickSort(std::vector<Sortable>& sortElements, int timeSleep, const std::atomic<bool>& interrupt);
	int heapSort(std::vector<Sortable>& sortElements, int timeSleep, const std::atomic<bool>& interrupt);
	int mergeSort(std::vector<Sortable>& sortElements, int timeSleep, const std::atomic<bool>& interrupt);
}

namespace algoUtils {
	void swap(std::vector<Sortable>& sortElements, int timeSleep, Sortable& el1, Sortable& el2);
}
