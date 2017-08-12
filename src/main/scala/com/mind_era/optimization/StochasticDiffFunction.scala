package com.mind_era.optimization

// Based on https://github.com/scalanlp/breeze/blob/releases/v0.13.1/math/src/main/scala/breeze/optimize/StochasticDiffFunction.scala

import spire.algebra.{Eq, Field}

/**
 * A differentiable function whose output is not guaranteed to be the same
 * across consecutive invocations.
 * @author dlwh
 */
trait StochasticDiffFunction[T, @specialized(Double) R] extends (T=>R)/* with NumericOps[StochasticDiffFunction[T]]*/ { outer =>
//  override def repr: StochasticDiffFunction[T] = this

  /** calculates the gradient at a point */
  def gradientAt(x: T): T = calculate(x)._2
  /** calculates the value at a point */
  def valueAt(x:T): R = calculate(x)._1

  final def apply(x:T): R = valueAt(x)

  /** Calculates both the value and the gradient at a point */
  def calculate(x:T):(R,T)

//  /**
//   * Lenses provide a way of mapping between two types, which we typically
//   * use to convert something to a DenseVector or other Tensor for optimization purposes.
//   */
//  def throughLens[U](implicit l: Isomorphism[T,U]):StochasticDiffFunction[U] = new StochasticDiffFunction[U] {
//    override def calculate(u: U) = {
//      val t = l.backward(u)
//      val (obj,gu) = outer.calculate(t)
//      (obj,l.forward(gu))
//    }
//  }
}