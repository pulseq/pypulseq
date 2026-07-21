# inserted for loop support by mveldmann

from types import SimpleNamespace
from typing import Union

from pypulseq.opts import Opts


def make_loop(
    loop_id: int,
    on_off: int = 0,
    system: Union[Opts, None] = None,
) -> SimpleNamespace:
    """
     Create a loop even

    See also `pypulseq.Sequence.sequence.Sequence.add_block()`.

    Parameters
    ----------
     loop_id : int
         ID of the loop.
     on_off : int, default=0
     system : Opts, default=Opts()
         System limits.

    Returns
    -------
     loop : SimpleNamespace
         loop event.

    Raises
    ------
     ValueError
         If invalid `channel` is passed. Must be one of 'physio1' or 'physio2'.
    """
    if system is None:
        system = Opts.default

    loop = SimpleNamespace()
    loop.type = 'loop'
    loop.loop_id = loop_id
    loop.on_off = on_off

    return loop
