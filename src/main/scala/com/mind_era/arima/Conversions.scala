package com.mind_era.arima

import spire.math.{Rational, Real}

import scala.languageFeature.implicitConversions

/**
  * Created by aborg on 09/06/2017.
  */
object Conversions {
  implicit class Id[A](val a: A) extends AnyVal with Convert[A, A] {
    override def to(a: A): A = a
  }

  implicit class RatToReal(val r: Rational) extends AnyVal with Convert[Rational, Real] {
    override def to(a: Rational): Real = a.toReal
  }

//  implicit class DoubleToRational(val d: Double) extends AnyVal with Convert[Double, Rational] {
//    override def to(a: Double): Rational = Rational(d)
//  }
//
//  implicit class DoubleToReal(val d: Double) extends AnyVal with Convert[Double, Real] {
//    override def to(a: Double): Real = Real(d)
//  }
//
  implicit def ratToReal: Convert[Rational, Real] = (a: Rational) => a.toReal

//  implicit def doubleToRational: Convert[Double, Rational] = (d: Double) => Rational(d)
//
//  implicit def doubleToReal: Convert[Double, Real] = (d: Double) => Real(d)

  implicit def orderToOrdering[T: spire.algebra.Order]: Ordering[T] = Ordering.fromLessThan(
    implicitly[spire.algebra.Order[T]].lt)

//  type ToDouble[A] = Convert[A, Double]
//  type FromDouble[B] = Convert[Double, B]
}

@deprecated("Use spire.math.ConvertableTo instead.")
trait Convert[A, B] extends Any {
  def to(a: A): B
}
