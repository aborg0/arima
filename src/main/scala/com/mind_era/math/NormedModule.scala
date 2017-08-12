package com.mind_era.math

//import spire.implicits._
import spire.algebra.{Field, Order, Signed, VectorSpace}
import spire.syntax.all._

// Based on https://github.com/scalanlp/breeze/blob/releases/v0.13.1/math/src/main/scala/breeze/math/VectorSpace.scala

//import com.mind_era.math.norm

trait Normed[V] {
  implicit def normImpl: norm.Impl[V, V]
}

trait NormedModule[V, S/*: Signed: Order: Field: Normed*/] extends VectorSpace[V, S] with Normed[S] {
  implicit def scalarSigned: Signed[S]
  implicit def scalarOrder: Order[S]
  implicit def scalarField: Field[S]
  implicit def scalarNormed: Normed[S]
  implicit def vectorNorm: norm.Impl[V, S]
  implicit def scalarNorm: norm.Impl[S, S] = norm.scalarNorm[S]
  def close(a: V, b: V, tolerance: S): Boolean = {
    norm[V, S](plus(a, negate(b))) <= tolerance *
      Order.max(scalar.one, Order.max(norm[V, S](a), norm[V, S](b)))
  }
}