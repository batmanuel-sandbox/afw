
// -*- LSST-C++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2018 LSST/AURA.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
#ifndef LSST_AFW_MATH_POLYNOMIALS_PackedInded_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_PackedIndex_h_INCLUDED

namespace lsst { namespace afw { namespace math { namespace polynomials {


/**
 *  A custom tuple that relates the indices of two 1-d functions for x and y
 *  to the flattened index for the 2-d function they form.
 *
 *  The packing algorithm is not specified by Index2d itself; it is intended
 *  to be a common data structure for other classes that define possibly
 *  different algorithms (e.g. PackedIndexIterator).
 */
struct Index2d {

    /// Construct an index with zero entrires.
    constexpr Index2d() noexcept : flat(0), x(0), y(0) {}

    /// Construct with the provided values.
    constexpr Index2d(std::size_t flat_, std::size_t x_, std::size_t y_) noexcept :
        flat(flat_), x(x_), y(y_)
    {}

    /// Equality comparison.
    constexpr bool operator==(Index2d const & other) const noexcept {
        return flat == other.flat && x == other.x && y == other.y;
    }

    /// Inequality comparison.
    constexpr bool operator!=(Index2d const & other) const noexcept {
        return !(*this == other);
    }

    std::size_t flat; ///< Index into the flattened 2-d function.
    std::size_t x;    ///< Index into the 1-d function for x.
    std::size_t y;    ///< Index into the 1-d functoin for y.
};

/**
 *  An iterator for traversing "packed" triangular 2-d series expansions,
 *  in which two 1-d expansions are truncated according to the sum of their
 *  orders and all values for one order are stored before any values of the
 *  subsequent order.
 *
 *  PackedIndexIterator dereferences to Index2d.  Typical usage is via a
 *  PackedIndexRange.
 *
 *  A pair of indices @f$(p, q)@f$ is mapped to the flattened position
 *  @f$i = (x+y)(x+y+1)/2 + x@f$, which yields the (x, y) ordering
 *  ```
 *      (0, 0), (0, 1), (1, 0), (1, 1), (0, 2), (1, 1), (2, 0), ...
 *  ```
 *  This packing ensures that the coefficients for an nth-order expansion are
 *  a contiguous subset of the coefficients for an (n+1)th-order expansion.
 *
 *  PackedIndexIterator is an STL input iterator, *not* a forward iterator;
 *  incrementing the iterator invalidates all previously-dereferenced values.
 */
class PackedIndexIterator {
public:

    using difference_type = std::ptrdiff_t;
    using value_type = Index2d;
    using pointer = Index2d const *;
    using reference = Index2d const &;
    using iterator_category = std::input_iterator_tag;

    /// Return the flattened offset to the start of the given order.
    static constexpr std::size_t computeOffset(std::size_t order) noexcept {
        return order*(order + 1)/2;
    }

    /// Return the flattened size of an expansion with the given maximum order (inclusive).
    static constexpr std::size_t computeSize(std::size_t order) noexcept {
        return computeOffset(order + 1);
    }

    /// Return the flattened index for the element with the given x and y orders.
    static constexpr std::size_t computeIndex(std::size_t x, std::size_t y) noexcept {
        return computeOffset(x + y) + x;
    }

    /// Construct an iterator one past the end of an expansion with the given order.
    static constexpr PackedIndexIterator makeEnd(std::size_t order) noexcept {
        return PackedIndexIterator(order);
    }

    /// Construct an iterator at the beginning of an expansion of any order.
    constexpr PackedIndexIterator() noexcept : _index() {}

    /// Construct an iterator pointing to the element with the given x and y orders.
    constexpr PackedIndexIterator(std::size_t x, std::size_t y) noexcept :
        _index(computeIndex(x, y), x, y)
    {}

    /// Dereference the iterator, yielding a Index2d const reference.
    constexpr reference operator*() const noexcept { return _index; }

    /// Dereference the iterator, yielding a Index2d const pointer.
    constexpr pointer operator->() const noexcept { return &_index; }

    /// Move to the next element in the packed array and return the iterator.
    PackedIndexIterator & operator++() noexcept {
        ++_index.flat;
        if (_index.y == 0) {
            _index.y = _index.x + 1;
            _index.x = 0;
        } else {
            --_index.y;
            ++_index.x;
        }
        return *this;
    }

    /// Move to the next element in the packed array and return a copy of the iterator before the move.
    PackedIndexIterator operator++(int) noexcept {
        PackedIndexIterator r(*this);
        ++(*this);
        return r;
    }

    /// Equality comparison.
    constexpr bool operator==(PackedIndexIterator const & other) const noexcept {
        return _index == other._index;
    }

    /// Inequality comparison.
    constexpr bool operator!=(PackedIndexIterator const & other) const noexcept {
        return !(*this == other);
    }

private:

    constexpr PackedIndexIterator(std::size_t order) noexcept :
        _index(computeOffset(order + 1), 0, order + 1)
    {}

    Index2d _index;
};

/**
 *  A specialized iterator range class for PackedIteratorRange, providing
 *  size calculation, comparison, and range-based `for` support.
 *
 *  @see PackedIndexIterator for information on the packing algorithm.
 */
class PackedIndexRange {
public:

    using value_type = PackedIndexIterator::value_type;
    using reference = PackedIndexIterator::reference;
    using pointer = PackedIndexIterator::pointer;
    using iterator = PackedIndexIterator;
    using const_iterator = PackedIndexIterator;
    using difference_type = PackedIndexIterator::difference_type;
    using size_type = std::size_t;

    /// Return the flattened offset to the start of the given order.
    static constexpr std::size_t computeOffset(std::size_t order) noexcept {
        return iterator::computeOffset(order);
    }

    /// Return the flattened size of an expansion with the given maximum order (inclusive).
    static constexpr std::size_t computeSize(std::size_t order) noexcept {
        return iterator::computeSize(order);
    }

    /// Return the flattened index for the element with the given x and y orders.
    static constexpr std::size_t computeIndex(std::size_t x, std::size_t y) noexcept {
        return iterator::computeIndex(x, y);
    }

    /// Construct from begin and end iterators.
    constexpr PackedIndexRange(iterator first, iterator last) noexcept :
        _begin(first),
        _end(last)
    {}

    /// Return an iterator to the start of the range.
    constexpr PackedIndexIterator begin() const noexcept { return _begin; }

    /// Return an iterator to the start of the range.
    constexpr PackedIndexIterator cbegin() const noexcept { return _begin; }

    /// Return an iterator to one past the end of the range.
    constexpr PackedIndexIterator end() const noexcept { return _end; }

    /// Return an iterator to one past the end of the range.
    constexpr PackedIndexIterator cend() const noexcept { return _end; }

    /// Return the number of elements in the flattened expansion.
    constexpr std::size_t size() const noexcept { return _end->flat - _begin->flat; }

    /// Return true if the number of elements in the flattened expansion is zero.
    constexpr bool empty() const noexcept { return size() == 0u; }

    /// Equality comparison.
    constexpr bool operator==(PackedIndexRange const & other) const noexcept {
        return _begin == other._begin && _end == other._end;
    }

    /// Inequality comparison.
    constexpr bool operator!=(PackedIndexRange const & other) const noexcept {
        return !(*this == other);
    }

private:
    PackedIndexIterator _begin;
    PackedIndexIterator _end;
};

}}}} // namespace lsst::afw::math::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_PackedIndex_h_INCLUDED
