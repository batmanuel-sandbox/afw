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
#ifndef LSST_AFW_MATH_POLYNOMIALS_Scaling1d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_Scaling1d_h_INCLUDED


namespace lsst { namespace afw { namespace math { namespace polynomials {

/**
 *  A 1-d affine transform that can be used to map one interval to another.
 *
 *  The transform is represented as an additive shift followed by a
 *  multiplicative scaling.
 *
 *  @note This class (and especially its 2-d counterpart, Scaling2d) has a lot
 *        in common with afw::geom::AffineTransform, and ideally they should
 *        share code and be highly interoperable.  Doing that well would
 *        require a larger-scale rethink of afw::geom, however, and at present
 *        we don't actually have a use case for that interoperability, so it's
 *        something we should keep in mind for the future, not a high priority
 *        for the present.
 *
 *  @exceptsafe All operations on Scaling1d are `noexcept`.
 *
 *  @see makeUnitRangeScaling1d
 */
class Scaling1d {
public:

    /// Construct from the given multiplicative scale and additive shift.
    Scaling1d(double scale, double shift) noexcept:
        _scale(scale),
        _shift(shift)
    {}

    /// Default copy constructor.
    Scaling1d(Scaling1d const &) = default;

    /// Default move constructor.
    Scaling1d(Scaling1d &&) = default;

    /// Default copy assignment.
    Scaling1d & operator=(Scaling1d const &) = default;

    /// Default move assignment.
    Scaling1d & operator=(Scaling1d &&) = default;

    /**
     *  Apply the transform in the forward direction.
     *
     *  Result is defined to be `(x + getShift()) * getScale()`.
     */
    double applyForward(double x) const noexcept {
        return (x + getShift())*getScale();
    }

    /**
     *  Apply the inverse of the forward transform;
     */
    double applyInverse(double y) const noexcept {
        return y/_scale - _shift;
    }

    /// Return the multiplicative scaling.
    double getScale() const noexcept { return _scale; }

    /// Return the additive shift.
    double getShift() const noexcept { return _shift; }

    /**
     *  Invert the transform.
     *
     *  If `r = t.invert()`, then `r.applyForward(x)` is equivalent to
     *  `t.applyInverse(x)` and `r.applyInverse(y)` is equivalent to
     *  `t.applyForward(y)`.
     */
    Scaling1d invert() const noexcept {
        return Scaling1d(1.0/_scale, -_shift*_scale);
    }

    /**
     *  Compose two transforms.
     *
     *  If `r = a.then(b)`, then `r.applyForward(x)` is equivalent to
     *  `b.applyForward(a.applyForward(x))` and `r.applyInverse(y)` is
     *  equivalent to `a.applyInverse(b.applyInverse(y))`.
     */
    Scaling1d then(Scaling1d const & second) const noexcept {
        return Scaling1d(getScale()*second.getScale(),
                                 getShift() + second.getShift()/getScale());
    }

private:
    double _scale;
    double _shift;
};

/**
 *  Return a Scaling1d that maps the interval [min, max] to [-1, 1].
 */
inline Scaling1d makeUnitRangeScaling1d(double min, double max) noexcept {
    return Scaling1d(2.0/(max - min), -0.5*(min + max));
}

}}}} // namespace lsst::afw::math::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_Scaling1d_h_INCLUDED
