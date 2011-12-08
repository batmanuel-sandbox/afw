// -*- lsst-c++ -*-
#ifndef AFW_TABLE_KeyBase_h_INCLUDED
#define AFW_TABLE_KeyBase_h_INCLUDED

namespace lsst { namespace afw { namespace table { 

template <typename T> class Key;
template <typename T> class Point;
template <typename T> class Shape;
template <typename T> class Array;
template <typename T> class Covariance;

/// @brief A base class for Key that allows subfield keys to be extracted.
template <typename T>
class KeyBase {};

/// @brief A base class for Key that allows subfield keys to be extracted.
template <typename U>
class KeyBase< Point<U> > {
public:
    Key<U> getX() const; ///< @brief Return a subfield key for the 'x' coordinate.
    Key<U> getY() const; ///< @brief Return a subfield key for the 'y' coordinate.
};

/// @brief A base class for Key that allows subfield keys to be extracted.
template <typename U>
class KeyBase< Shape<U> > {
public:
    Key<U> getIXX() const; ///< @brief Return a subfield key for the 'xx' value.
    Key<U> getIYY() const; ///< @brief Return a subfield key for the 'yy' value.
    Key<U> getIXY() const; ///< @brief Return a subfield key for the 'xy' value.
};

/// @brief A base class for Key that allows subfield keys to be extracted.
template <typename U>
class KeyBase< Array<U> > {
public:
    Key<U> operator[](int i) const; ///< @brief Return a subfield key for the i-th element of the array.
};

/// @brief A base class for Key that allows subfield keys to be extracted.
template <typename U>
class KeyBase< Covariance<U> > {
public:
    ///< @brief Return a subfield key for element (i,j) of the covariance matrix.
    Key<U> operator()(int i, int j) const;
};

/// @brief A base class for Key that allows subfield keys to be extracted.
template <typename U>
class KeyBase< Covariance< Point<U> > > {
public:
    ///< @brief Return a subfield key for element (i,j) of the covariance matrix.
    Key<U> operator()(int i, int j) const;
};

/// @brief A base class for Key that allows subfield keys to be extracted.
template <typename U>
class KeyBase< Covariance< Shape<U> > > {
public:
    ///< @brief Return a subfield key for element (i,j) of the covariance matrix.
    Key<U> operator()(int i, int j) const;
};

}}} // namespace lsst::afw::table

#endif // !AFW_TABLE_KeyBase_h_INCLUDED
