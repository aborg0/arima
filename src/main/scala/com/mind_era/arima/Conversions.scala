package com.mind_era.arima

import spire.math.{Rational, Real}

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

  implicit def ratToReal: Convert[Rational, Real] = (a: Rational) => a.toReal
}

trait Convert[A, B] extends Any {
  def to(a: A): B
}
