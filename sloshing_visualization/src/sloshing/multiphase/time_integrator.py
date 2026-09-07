"""Explicit BE startup schedules and mathematically fixed-spacing main BDF2."""
import math


def derivative_coefficients(step, dt, scheme="bdf2"):
    if step < 1 or dt <= 0 or scheme not in ("be", "bdf2"):
        raise ValueError("Invalid time integration arguments")
    if step == 1 or scheme == "be":
        return 1 / dt, -1 / dt, 0.0
    return 1.5 / dt, -2 / dt, 0.5 / dt
class IntegrationSchedule:
    """Explicit BE startup blocks, one main BE history interval, then fixed BDF2.

    No BDF2 derivative ever uses unequal spacing. Block counts are integral and
    checked before assembly; accepted endpoint times are computed, not accumulated.
    """
    def __init__(self,config):
        self.blocks=[]
        start=0.;offset=0
        startup=list(config.startup_stages)
        if config.startup_dt is not None:
            startup=[{"dt":config.startup_dt,"t_end":config.startup_t_end}]
        for block in startup:
            count=self._count(block["t_end"]-start,block["dt"])
            self.blocks.append({"start":start,"end":block["t_end"],"dt":block["dt"],
                "count":count,"offset":offset,"phase":"startup_be","scheme":"be"})
            start=block["t_end"];offset+=count
        count=self._count(config.t_end-start,config.dt)
        if count:
            self.blocks.append({"start":start,"end":config.t_end,"dt":config.dt,
                "count":count,"offset":offset,"phase":"main","scheme":config.time_scheme})
        self.nsteps=offset+count
        self.descriptor={"blocks":self.blocks,"total_steps":self.nsteps,
            "switch_rule":"first main interval BE, then fixed-spacing BDF2; no micro-history reuse"}

    @staticmethod
    def _count(duration,dt):
        count=round(duration/dt)
        if count<0 or not math.isclose(count*dt,duration,rel_tol=1e-10,abs_tol=1e-14):
            raise ValueError("Each schedule block must contain an integral number of its dt intervals")
        return count

    def step(self,accepted_steps):
        if not 0<=accepted_steps<self.nsteps:
            raise ValueError("No further interval in the declared schedule")
        for block in self.blocks:
            local=accepted_steps-block["offset"]
            if 0<=local<block["count"]:
                if block["phase"]=="startup_be":
                    phase,scheme="startup_be","be"
                elif block["scheme"]=="be":
                    phase,scheme="main_be","be"
                elif local==0:
                    phase,scheme="main_be_history","be"
                else:
                    phase,scheme="bdf2","bdf2"
                return {"dt":block["dt"],"time":block["start"]+(local+1)*block["dt"],
                        "phase":phase,"scheme":scheme}
        raise RuntimeError("Schedule block not found")
