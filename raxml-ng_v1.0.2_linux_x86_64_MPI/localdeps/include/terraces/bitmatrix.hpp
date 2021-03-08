#ifndef TERRACES_BITMATRIX_HPP
#define TERRACES_BITMATRIX_HPP

#include <vector>

#include "trees.hpp"

namespace terraces {

/** A (memory-wise) compact bitmatrix. */
class bitmatrix {
public:
	/** Constructs a bitmatrix with \p rows rows and \p cols columns. */
	bitmatrix(index_t rows, index_t cols);

	/** @returns the number of rows. */
	index_t rows() const;
	/** @returns the number of columns. */
	index_t cols() const;

	/** Returns the entry at cordinates (row, col). */
	bool get(index_t row, index_t col) const;
	/** Sets the entry at coordinates (row, col). */
	void set(index_t row, index_t col, bool val);
	/** Writes the bit rows \p in1 | \p in2 to row \p out. */
	void row_or(index_t in1, index_t in2, index_t out);

	/** Returns a bitmatrix containing only the given columns. */
	bitmatrix get_cols(const std::vector<std::size_t>& cols) const;

	/** Tests for equality with another bitmatrix. */
	bool operator==(const bitmatrix& other) const;
	/** Tests for inequality with another bitmatrix. */
	bool operator!=(const bitmatrix& other) const;

private:
	index_t m_rows;
	index_t m_cols;
	std::vector<bool> m_vec; // yes, it's evil
};
} // namespace terraces

#endif // TERRACES_BITMATRIX_HPP
