package com.mind_era.math

import spire.algebra.{Eq, Field, Signed}

// Based on https://github.com/scalanlp/breeze/blob/releases/v0.13.1/math/src/main/scala/breeze/linalg/functions/norm.scala

//import breeze.generic.UFunc
//import breeze.macros.expand

/**
 * Computes the norm of an object. Many tensor objects have a norm implementation implicit, which is what this calls.
 */
object norm extends UFunc {
//  @expand
//  @expand.valify
  implicit def scalarNorm[/*@expand.args(Int, Long, Float, Double) */T: Signed: Field]: Impl[T, T] = new Impl[T, T] {
    def apply(v1: T): T = Signed.apply.abs(v1)
  }

  implicit def normalNormToNormUnit[T: Signed: Field](implicit normImpl: Impl[T, T]):Impl2[T, Unit, T] = {
    new Impl2[T, Unit, T] {
      def apply(v: T, x: Unit): T = normImpl(v)
    }
  }



//  implicit def normDoubleToNormalNorm[T: Signed: Field](implicit normImpl: Impl2[T, T, T]):Impl[T, T] = {
//    new Impl[T, T] {
//      def apply(v: T): T = normImpl(v, 2.0)
//    }
//  }
//
//  implicit def fromCanNormInt[T: Signed: Field](implicit impl: Impl2[T, T, T]):Impl2[T, Int, T] = {
//    new Impl2[T, Int, T] {
//      def apply(v: T, v2: Int): T = impl(v, v2)
//    }
//  }

  implicit def fromCanNormFloat[T](implicit impl: Impl2[T, Double, Double]):Impl2[T, Float, Float] = {
    new Impl2[T, Float, Float] {
      def apply(v: T, v2: Float): Float = impl(v, v2).toFloat
    }
  }
}