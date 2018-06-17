package com.mind_era.arima

/**
  * Created by aborg on 04/06/2017.
  */
final case class SingleValue[V] private(private val v: V) extends IndexedSeq[V] {
  override def length: Int = 1

  override def apply(idx: Int): V = {
    require(idx == 0, s"Only single value is present, but requested: $idx")
    v
  }
}
