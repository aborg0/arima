package com.mind_era

package object util extends IteratorImplicits {
}

trait IteratorImplicits {

  class RichIterator[T](iter: Iterator[T]) {
    def tee(f: T => Unit): Iterator[T] = new Iterator[T] {
      def next = {
        val n = iter.next;
        f(n);
        n
      }

      def hasNext = {
        iter.hasNext;
      }
    }

    def takeUpToWhere(f: T => Boolean): Iterator[T] = new Iterator[T] {
      var done = false

      def next = {
        if (done) throw new NoSuchElementException()
        val n = iter.next;
        done = f(n)
        n
      }

      def hasNext = {
        !done && iter.hasNext;
      }
    }

    def last = {
      var x = iter.next()
      while (iter.hasNext) {
        x = iter.next()
      }
      x
    }
  }

  implicit def scEnrichIterator[T](iter: Iterator[T]) = new RichIterator[T](iter)
}
