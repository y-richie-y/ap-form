# AP form simulator

Very basic Clifford state simulator using AP form in ~ 100 loc,
written for (my) educational purposes.

<img width="553" alt="image" src="https://github.com/user-attachments/assets/6412bddd-3bc1-479d-ba6b-7bb46033c6e7" />

Clifford states can be expressed as a list of projectors onto an affine subspace $A$ over $\mathbb{F}^n_2$, transformed by a quadratic phase polynomial $\phi(\vec{x})$.

You can read more about them in the [Picturing Quantum Software](https://github.com/zxcalc/book) book.

## Example
```python
from state import APState

APState(3).s(0).h(1).cx(0, 1).cz(0, 2).vec()
# APState(n_qubits=3, projectors=set(), phases4=[3, 3, 3], czs={CZ({1, 2}), CZ({0, 1}), CZ({0, 2})})
# array([ 0.5+0.j ,  0. +0.j ,  0. +0.j ,  0. +0.5j,  0.5+0.j ,  0. +0.j ,
#        0. +0.j , -0. -0.5j])
```
