**Added:**

* Add ``make_equation``, ``register_constant``, ``register_argument``,
  ``register_operator``, ``register_builder`` and ``deregister_builder`` to
  ``EquationFactory``, and ``wrap_argument``, ``wrap_operator``,
  ``wrap_function`` and ``get_builder`` to ``diffpy.srfit.equation.builder``.
* Add ``Equation.set_root``, ``Operator.add_literal`` and ``make_operator`` to
  ``diffpy.srfit.equation``.
* Add ``get_args``, ``get_expression`` and ``pretty_print`` to
  ``diffpy.srfit.equation.visitors``, and ``on_argument``, ``on_operator`` and
  ``on_equation`` to ``Visitor``.

**Changed:**

* ``OperatorABC`` now requires ``add_literal`` instead of ``addLiteral``, and
  ``Visitor`` subclasses now override ``on_argument``, ``on_operator`` and
  ``on_equation`` instead of ``onArgument``, ``onOperator`` and
  ``onEquation``.

**Deprecated:**

* ``EquationFactory.makeEquation``, ``registerConstant``, ``registerArgument``,
  ``registerOperator``, ``registerBuilder`` and ``deRegisterBuilder``,
  ``wrapArgument``, ``wrapOperator``, ``wrapFunction``, ``getBuilder``,
  ``Equation.setRoot``, ``Operator.addLiteral``, ``makeOperator``,
  ``getArgs``, ``getExpression``, ``prettyPrint`` and ``Visitor.onArgument``,
  ``onOperator`` and ``onEquation`` are deprecated in favor of their
  snake_case replacements. The old names still work and warn, and will be
  removed in version 4.0.0.

**Removed:**

* <news item>

**Fixed:**

* Fix the test check for observers left on global equation builders.

**Security:**

* <news item>
