package com.mind_era.math

import spire.algebra.Order
import spire.math.{Rational, Real}

abstract class Bounds[V: Order] {
  def lowerBound: V
  def upperBound: V
}

object Bounds {
  implicit def DoubleBounds(implicit ev: Order[Double]): Bounds[Double] = new Bounds[Double] {
    override val upperBound: Double = ev.max(Double.PositiveInfinity, Double.NegativeInfinity)

    override val lowerBound: Double = ev.min(Double.PositiveInfinity, Double.NegativeInfinity)
  }

  /** Just approximation, can hold more extreme values too. */
  implicit def RealBounds(implicit ev: Order[Real]): Bounds[Real] = new Bounds[Real] {
    override val lowerBound: Real = ev.min(Real(Double.NegativeInfinity), Real(Double.PositiveInfinity))

    override val upperBound: Real = ev.max(Real(Double.NegativeInfinity), Real(Double.PositiveInfinity))
  }
  /** Just approximation, can hold more extreme values too. */
  implicit def RationalBounds(implicit ev: Order[Rational]): Bounds[Rational] = new Bounds[Rational] {
    override val lowerBound: Rational = ev.min(Rational(Double.NegativeInfinity), Rational(Double.PositiveInfinity))

    override val upperBound: Rational = ev.max(Rational(Double.NegativeInfinity), Rational(Double.PositiveInfinity))
  }
}


