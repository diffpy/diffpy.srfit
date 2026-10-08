"""Unit tests for observable.py."""

from diffpy.srfit.util.observable import Observable


class _Recorder:
    """Observer that records the arguments it is notified with."""

    def __init__(self):
        self.calls = []

    def record(self, semaphores):
        self.calls.append(semaphores)


def test_add_and_remove_observer():
    # C1: An observer is added, then removed.
    # Expected: has_observer is True after adding and False after removing.
    observable = Observable()
    recorder = _Recorder()
    assert observable.has_observer(recorder.record) is False
    observable.add_observer(recorder.record)
    assert observable.has_observer(recorder.record) is True
    observable.remove_observer(recorder.record)
    actual_has_observer = observable.has_observer(recorder.record)
    expected_has_observer = False
    assert actual_has_observer == expected_has_observer


def test_notify():
    # C1: Two observers are added and the observable notifies with an extra
    # object.
    # Expected: Each observer is called once with the observable followed
    # by the extra object.
    observable = Observable()
    first, second = _Recorder(), _Recorder()
    observable.add_observer(first.record)
    observable.add_observer(second.record)
    observable.notify(("extra",))
    expected_calls = [(observable, "extra")]
    assert first.calls == expected_calls
    assert second.calls == expected_calls
