#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <atomic>
/*
 * Multiplication in GF(2^8)
 * http://en.wikipedia.org/wiki/Finite_field_arithmetic
 * Irreducible polynomial m(x) = x8 + x4 + x3 + x + 1
 * 
 * NOTE: we are using the look up table instead of the (slower) gmult function
 */
#define gmult(a,b) gmult_aes[256*a + b]

extern uint8_t gmult_aes[];

uint8_t *aes_init(size_t key_size);

void aes_key_expansion(uint8_t *key, uint8_t *w);

void aes_inv_cipher(uint8_t *in, uint8_t *out, uint8_t *w);

void aes_cipher(uint8_t *in, uint8_t *out, uint8_t *w);

class SpinLock {
 
public:
 SpinLock() : flag_(false)
 {}
 
 bool lock()
 {
 bool expect = false;
 while (!flag_.compare_exchange_weak(expect, true))
 {
  //这里一定要将expect复原，执行失败时expect结果是未定的
  expect = false;
  return false;
 }
 return true;
 }
 
 void unlock()
 {
 flag_.store(false);
 }
 
private:
 std::atomic<bool> flag_;
};