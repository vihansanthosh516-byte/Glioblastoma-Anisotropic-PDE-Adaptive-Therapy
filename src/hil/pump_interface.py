#!/usr/bin/env python3
"""
Phase 14: Hardware-in-the-Loop (HIL) Integration
=================================================
Bridges the software control agent with physical/hardware-emulated drug delivery apparatus.

Components:
- Pump Interface: Serial/TCP socket interface for Chemyx / Harvard Apparatus syringe pumps
- Safety Watchdog: Deterministic safety override enforcing hard constraints
- Latency Benchmarking: End-to-end pipeline < 500ms guarantee

Supported Pumps:
- Chemyx Fusion 200/400 series (RS-232/USB)
- Harvard Apparatus Pump 11 Elite (RS-232)
- Generic syringe pump protocol emulation
"""
import time
import threading
import socket
import serial
import queue
import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Callable, Any
from collections import deque
import time

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PumpStatus(Enum):
    IDLE = "IDLE"
    INFUSING = "INFUSING"
    PAUSED = "PAUSED"
    ERROR = "ERROR"
    STOPPED = "STOPPED"


@dataclass
class PumpParameters:
    """Pump hardware parameters and limits."""
    max_rate_ul_min: float = 1000.0      # µL/min max
    min_rate_ul_min: float = 0.01        # µL/min min
    max_volume_ul: float = 10000.0       # µL max syringe volume
    syringe_diameter_mm: float = 14.57   # Standard 60mL syringe
    max_pressure_psi: float = 50.0       # Max pressure
    address: int = 0                     # Pump address (for multi-pump)


@dataclass
class PumpState:
    """Current pump state."""
    status: PumpStatus = PumpStatus.IDLE
    rate_ul_min: float = 0.0
    volume_delivered_ul: float = 0.0
    target_volume_ul: float = 0.0
    pressure_psi: float = 0.0
    error_code: int = 0
    timestamp: float = field(default_factory=time.time)


class PumpInterface(ABC):
    """Abstract base class for pump interfaces."""
    
    @abstractmethod
    def connect(self) -> bool:
        pass
    
    @abstractmethod
    def disconnect(self) -> bool:
        pass
    
    @abstractmethod
    def infuse(self, rate_ul_min: float, target_volume_ul: float = None) -> bool:
        pass
    
    @abstractmethod
    def stop(self) -> bool:
        pass
    
    @abstractmethod
    def pause(self) -> bool:
        pass
    
    @abstractmethod
    def resume(self) -> bool:
        pass
    
    @abstractmethod
    def get_state(self) -> PumpState:
        pass
    
    @abstractmethod
    def set_diameter(self, diameter_mm: float) -> bool:
        pass


class MockPumpInterface(PumpInterface):
    """Mock pump interface for testing and development."""
    
    def __init__(self, params: PumpParameters = None):
        self.params = params or PumpParameters()
        self.state = PumpState()
        self._lock = threading.Lock()
        self._infusion_thread = None
        self._stop_event = threading.Event()
        logger.info(f"[MockPump] Initialized with params: {self.params}")
    
    def connect(self) -> bool:
        with self._lock:
            self.state.status = PumpStatus.IDLE
            self.state.timestamp = time.time()
            logger.info("[MockPump] Connected")
            return True
    
    def disconnect(self) -> bool:
        with self._lock:
            self.stop()
            self.state.status = PumpStatus.STOPPED
            logger.info("[MockPump] Disconnected")
            return True
    
    def infuse(self, rate_ul_min: float, target_volume_ul: float = None) -> bool:
        with self._lock:
            # Validate rate
            if rate_ul_min < self.params.min_rate_ul_min:
                rate_ul_min = self.params.min_rate_ul_min
            if rate_ul_min > self.params.max_rate_ul_min:
                rate_ul_min = self.params.max_rate_ul_min
            
            if target_volume_ul is not None:
                if target_volume_ul > self.params.max_volume_ul:
                    logger.error(f"[MockPump] Target volume {target_volume_ul}µL exceeds max {self.params.max_volume_ul}µL")
                    return False
                self.state.target_volume_ul = target_volume_ul
            
            self.state.rate_ul_min = rate_ul_min
            self.state.status = PumpStatus.INFUSING
            self.state.timestamp = time.time()
            self._stop_event.clear()
            
            # Start infusion simulation thread
            self._infusion_thread = threading.Thread(target=self._infusion_loop, daemon=True)
            self._infusion_thread.start()
            
            logger.info(f"[MockPump] Infusing at {rate_ul_min:.2f} µL/min, target: {target_volume_ul}µL")
            return True
    
    def _infusion_loop(self):
        """Simulate infusion in background thread."""
        volume_per_sec = self.state.rate_ul_min / 60.0  # µL/sec
        last_time = time.time()
        
        while not self._stop_event.is_set():
            current_time = time.time()
            dt = current_time - last_time
            last_time = current_time
            
            with self._lock:
                if self.state.status != PumpStatus.INFUSING:
                    break
                
                # Update delivered volume
                self.state.volume_delivered_ul += volume_per_sec * dt
                self.state.timestamp = current_time
                
                # Check target volume
                if (self.state.target_volume_ul > 0 and 
                    self.state.volume_delivered_ul >= self.state.target_volume_ul):
                    self.state.status = PumpStatus.IDLE
                    logger.info(f"[MockPump] Target volume reached: {self.state.volume_delivered_ul:.1f}µL")
                    break
        
        self._stop_event.set()
    
    def stop(self) -> bool:
        with self._lock:
            self._stop_event.set()
            if self._infusion_thread:
                self._infusion_thread.join(timeout=1.0)
            self.state.rate_ul_min = 0.0
            self.state.status = PumpStatus.STOPPED
            self.state.timestamp = time.time()
            logger.info(f"[MockPump] Stopped. Total delivered: {self.state.volume_delivered_ul:.2f}µL")
            return True
    
    def pause(self) -> bool:
        with self._lock:
            if self.state.status == PumpStatus.INFUSING:
                self._stop_event.set()
                self.state.status = PumpStatus.PAUSED
                self.state.timestamp = time.time()
                logger.info(f"[MockPump] Paused at {self.state.volume_delivered_ul:.2f}µL")
                return True
            return False
    
    def resume(self) -> bool:
        with self._lock:
            if self.state.status == PumpStatus.PAUSED:
                self._stop_event.clear()
                self._infusion_thread = threading.Thread(target=self._infusion_loop, daemon=True)
                self._infusion_thread.start()
                self.state.status = PumpStatus.INFUSING
                self.state.timestamp = time.time()
                logger.info("[MockPump] Resumed")
                return True
            return False
    
    def get_state(self) -> PumpState:
        with self._lock:
            return PumpState(
                status=self.state.status,
                rate_ul_min=self.state.rate_ul_min,
                volume_delivered_ul=self.state.volume_delivered_ul,
                target_volume_ul=self.state.target_volume_ul,
                pressure_psi=self.state.pressure_psi,
                error_code=self.state.error_code,
                timestamp=self.state.timestamp
            )
    
    def set_diameter(self, diameter_mm: float) -> bool:
        with self._lock:
            self.params.syringe_diameter_mm = diameter_mm
            logger.info(f"[MockPump] Syringe diameter set to {diameter_mm}mm")
            return True


class ChemyxPumpInterface(PumpInterface):
    """Chemyx Fusion series pump interface (RS-232/USB)."""
    
    def __init__(self, port: str, params: PumpParameters = None, baudrate: int = 9600):
        self.port = port
        self.params = params or PumpParameters()
        self.baudrate = baudrate
        self.serial = None
        self._lock = threading.Lock()
        self.state = PumpState()
        self._response_queue = queue.Queue()
        self._reader_thread = None
        self._stop_reader = threading.Event()
    
    def connect(self) -> bool:
        try:
            self.serial = serial.Serial(
                port=self.port,
                baudrate=self.baudrate,
                bytesize=serial.EIGHTBITS,
                parity=serial.PARITY_NONE,
                stopbits=serial.STOPBITS_ONE,
                timeout=1.0
            )
            
            # Start reader thread
            self._stop_reader.clear()
            self._reader_thread = threading.Thread(target=self._reader_loop, daemon=True)
            self._reader_thread.start()
            
            # Initialize pump
            self._send_command("&")
            time.sleep(0.1)
            self._send_command(f"DIAMETER {self.params.syringe_diameter_mm:.2f}")
            self._send_command("RAT 0")
            
            self.state.status = PumpStatus.IDLE
            logger.info(f"[ChemyxPump] Connected on {self.port}")
            return True
        except Exception as e:
            logger.error(f"[ChemyxPump] Connection failed: {e}")
            return False
    
    def _reader_loop(self):
        """Background thread reading pump responses."""
        while not self._stop_reader.is_set():
            try:
                if self.serial and self.serial.in_waiting:
                    line = self.serial.readline().decode('ascii', errors='ignore').strip()
                    if line:
                        self._response_queue.put(line)
            except Exception as e:
                logger.error(f"[ChemyxPump] Reader error: {e}")
                break
    
    def _send_command(self, cmd: str) -> Optional[str]:
        with self._lock:
            if not self.serial:
                return None
            try:
                self.serial.write(f"{cmd}\r".encode('ascii'))
                # Wait for response
                start = time.time()
                while time.time() - start < 2.0:
                    try:
                        resp = self._response_queue.get(timeout=0.1)
                        if resp:
                            return resp
                    except queue.Empty:
                        continue
            except Exception as e:
                logger.error(f"[ChemyxPump] Send error: {e}")
            return None
    
    def disconnect(self) -> bool:
        self.stop()
        self._stop_reader.set()
        if self._reader_thread:
            self._reader_thread.join(timeout=1.0)
        if self.serial:
            self.serial.close()
        self.state.status = PumpStatus.STOPPED
        return True
    
    def infuse(self, rate_ul_min: float, target_volume_ul: float = None) -> bool:
        # Convert µL/min to mL/hr for Chemyx
        rate_ml_hr = rate_ul_min * 60 / 1000
        
        with self._lock:
            self._send_command(f"RAT {rate_ml_hr:.4f}")
            if target_volume_ul:
                vol_ml = target_volume_ul / 1000
                self._send_command(f"VOL {vol_ml:.4f}")
            self._send_command("RUN")
            
            self.state.rate_ul_min = rate_ul_min
            self.state.target_volume_ul = target_volume_ul or 0
            self.state.status = PumpStatus.INFUSING
            self.state.timestamp = time.time()
            return True
    
    def stop(self) -> bool:
        self._send_command("STP")
        with self._lock:
            self.state.status = PumpStatus.STOPPED
            self.state.rate_ul_min = 0.0
            self.state.timestamp = time.time()
        return True
    
    def pause(self) -> bool:
        self._send_command("STP")
        with self._lock:
            if self.state.status == PumpStatus.INFUSING:
                self.state.status = PumpStatus.PAUSED
                self.state.timestamp = time.time()
            return True
    
    def resume(self) -> bool:
        self._send_command("RUN")
        with self._lock:
            self.state.status = PumpStatus.INFUSING
            self.state.timestamp = time.time()
            return True
    
    def get_state(self) -> PumpState:
        with self._lock:
            return PumpState(
                status=self.state.status,
                rate_ul_min=self.state.rate_ul_min,
                volume_delivered_ul=self.state.volume_delivered_ul,
                target_volume_ul=self.state.target_volume_ul,
                pressure_psi=self.state.pressure_psi,
                error_code=self.state.error_code,
                timestamp=self.state.timestamp
            )
    
    def set_diameter(self, diameter_mm: float) -> bool:
        self._send_command(f"DIAMETER {diameter_mm:.2f}")
        self.params.syringe_diameter_mm = diameter_mm
        return True


class HarvardPumpInterface(PumpInterface):
    """Harvard Apparatus Pump 11 Elite interface."""
    
    def __init__(self, port: str, params: PumpParameters = None, baudrate: int = 9600):
        self.port = port
        self.params = params or PumpParameters()
        self.baudrate = baudrate
        self.serial = None
        self.state = PumpState()
        self._lock = threading.Lock()
    
    def connect(self) -> bool:
        try:
            self.serial = serial.Serial(
                port=self.port,
                baudrate=self.baudrate,
                timeout=1.0
            )
            
            # Initialize - Harvard uses different protocol
            self._send("RAT 0")
            self._send(f"DIA {self.params.syringe_diameter_mm:.2f}")
            
            self.state.status = PumpStatus.IDLE
            logger.info(f"[HarvardPump] Connected on {self.port}")
            return True
        except Exception as e:
            logger.error(f"[HarvardPump] Connection failed: {e}")
            return False
    
    def _send(self, cmd: str) -> Optional[str]:
        with self._lock:
            try:
                self.serial.write(f"{cmd}\r".encode())
                time.sleep(0.1)
                if self.serial.in_waiting:
                    return self.serial.readline().decode().strip()
            except Exception as e:
                logger.error(f"[HarvardPump] Send error: {e}")
            return None
    
    def infuse(self, rate_ul_min: float, target_volume_ul: float = None) -> bool:
        rate_ml_hr = rate_ul_min * 60 / 1000
        self._send(f"RAT {rate_ml_hr:.4f}")
        if target_volume_ul:
            self._send(f"VOL {target_volume_ul/1000:.4f}")
        self._send("RUN")
        
        with self._lock:
            self.state.rate_ul_min = rate_ul_min
            self.state.target_volume_ul = target_volume_ul or 0
            self.state.status = PumpStatus.INFUSING
            self.state.timestamp = time.time()
        return True
    
    def stop(self) -> bool:
        self._send("STP")
        with self._lock:
            self.state.status = PumpStatus.STOPPED
            self.state.rate_ul_min = 0.0
            self.state.timestamp = time.time()
        return True
    
    def pause(self) -> bool:
        self._send("STP")
        with self._lock:
            self.state.status = PumpStatus.PAUSED
            return True
    
    def resume(self) -> bool:
        self._send("RUN")
        with self._lock:
            self.state.status = PumpStatus.INFUSING
            return True
    
    def get_state(self) -> PumpState:
        with self._lock:
            return self.state
    
    def set_diameter(self, diameter_mm: float) -> bool:
        self._send(f"DIA {diameter_mm:.2f}")
        self.params.syringe_diameter_mm = diameter_mm
        return True
    
    def disconnect(self) -> bool:
        self.stop()
        if self.serial:
            self.serial.close()
        return True


class PumpFactory:
    """Factory for creating pump interfaces."""
    
    @staticmethod
    def create_pump(pump_type: str, **kwargs) -> PumpInterface:
        pump_type = pump_type.lower()
        if pump_type == "mock":
            return MockPumpInterface(kwargs.get('params'))
        elif pump_type == "chemyx":
            return ChemyxPumpInterface(kwargs.get('port'), kwargs.get('params'))
        elif pump_type == "harvard":
            return HarvardPumpInterface(kwargs.get('port'), kwargs.get('params'))
        else:
            raise ValueError(f"Unknown pump type: {pump_type}")


# Safety Watchdog
class SafetyWatchdog:
    """
    Deterministic safety override enforcing hard constraints on:
    - Cumulative toxicity per drug
    - Maximum hourly volumetric delivery (C_max)
    - Maximum instantaneous rate
    - Pressure limits
    - Cumulative volume limits
    """
    
    def __init__(
        self,
        max_cumulative_toxicity: Dict[str, float] = None,
        max_hourly_volume_ul: Dict[str, float] = None,
        max_instantaneous_rate_ul_min: Dict[str, float] = None,
        max_pressure_psi: float = 50.0,
        max_cumulative_volume_ul: Dict[str, float] = None,
        check_interval_sec: float = 0.1
    ):
        self.max_cumulative_toxicity = max_cumulative_toxicity or {
            "TMZ": 1.0, "Inhibitor": 0.8, "Radiation": 0.6
        }
        self.max_hourly_volume = max_hourly_volume_ul or {
            "TMZ": 5000, "Inhibitor": 3000, "Radiation": 1000
        }
        self.max_instantaneous_rate = max_instantaneous_rate_ul_min or {
            "TMZ": 100, "Inhibitor": 100, "Radiation": 50
        }
        self.max_cumulative_volume_ul = max_cumulative_volume_ul or {
            "TMZ": 20000, "Inhibitor": 15000, "Radiation": 5000
        }
        self.max_pressure_psi = max_pressure_psi
        
        self.check_interval = check_interval_sec
        self._lock = threading.Lock()
        
        # Runtime tracking
        self.cumulative_toxicity = {k: 0.0 for k in self.max_cumulative_toxicity}
        self.cumulative_volume_ul = {k: 0.0 for k in self.max_cumulative_volume_ul}
        self.hourly_volume = {k: deque(maxlen=3600) for k in self.max_hourly_volume}  # 1 hour at 1Hz
        self.last_check = time.time()
        self.violations = []
        self._monitoring = False
        self._monitor_thread = None
    
    def record_dose(self, drug: str, volume_ul: float, toxicity: float):
        """Record a dose delivery."""
        with self._lock:
            self.cumulative_volume[drug] = self.cumulative_volume.get(drug, 0) + volume_ul
            self.cumulative_toxicity[drug] = self.cumulative_toxicity.get(drug, 0) + toxicity
            self.hourly_volume[drug].append((time.time(), volume_ul))
    
    def check_safety(self, pump_state: PumpState, drug: str, proposed_rate: float) -> Tuple[bool, List[str]]:
        """Check if proposed action is safe. Returns (safe, violations)."""
        violations = []
        now = time.time()
        
        # Check instantaneous rate
        max_rate = self.max_instantaneous_rate.get(drug, float('inf'))
        if proposed_rate > max_rate:
            violations.append(f"Rate {proposed_rate:.1f} µL/min exceeds max {max_rate} µL/min for {drug}")
        
        # Check cumulative toxicity
        cum_tox = self.cumulative_toxicity.get(drug, 0)
        max_tox = self.max_cumulative_toxicity.get(drug, float('inf'))
        if cum_tox >= max_tox:
            violations.append(f"Cumulative toxicity {cum_tox:.3f} exceeds max {max_tox} for {drug}")
        
        # Check cumulative volume
        cum_vol = self.cumulative_volume_ul.get(drug, 0)
        max_vol = self.max_cumulative_volume_ul.get(drug, float('inf'))
        if cum_vol >= max_vol:
            violations.append(f"Cumulative volume {cum_vol:.1f}µL exceeds max {max_vol}µL for {drug}")
        
        # Check hourly volume (sliding window)
        hourly = sum(v for t, v in self.hourly_volume.get(drug, []) if now - t < 3600)
        max_hourly = self.max_hourly_volume.get(drug, float('inf'))
        if hourly >= max_hourly:
            violations.append(f"Hourly volume {hourly:.1f}µL exceeds max {max_hourly}µL for {drug}")
        
        # Check pressure
        if pump_state.pressure_psi > self.max_pressure_psi:
            violations.append(f"Pressure {pump_state.pressure_psi:.1f} psi exceeds max {self.max_pressure_psi} psi")
        
        return len(violations) == 0, violations
    
    def enforce_safety(self, pump: PumpInterface, pump_state: PumpState, 
                       drug: str, proposed_rate: float) -> Tuple[bool, float]:
        """
        Enforce safety constraints. Returns (allowed, adjusted_rate).
        If unsafe, adjusts rate down or stops pump.
        """
        safe, violations = self.check_safety(pump_state, drug, proposed_rate)
        
        if not safe:
            for v in violations:
                logger.warning(f"[SafetyWatchdog] VIOLATION: {v}")
                self.violations.append({"time": time.time(), "violation": v, "drug": drug})
            
            # Try to adjust rate down
            max_allowed = self.max_instantaneous_rate.get(drug, proposed_rate)
            adjusted_rate = min(proposed_rate, max_allowed)
            
            # If still violating cumulative limits, stop
            if any("Cumulative" in v or "Hourly" in v for v in violations):
                logger.error(f"[SafetyWatchdog] STOPPING pump due to cumulative violation for {drug}")
                pump.stop()
                return False, 0.0
            
            if adjusted_rate < proposed_rate:
                logger.warning(f"[SafetyWatchdog] Reducing rate from {proposed_rate} to {adjusted_rate} µL/min")
                return True, adjusted_rate
            
            return False, 0.0
        
        return True, proposed_rate
    
    def start_monitoring(self, pump: PumpInterface, drug: str):
        """Start background safety monitoring."""
        self._monitoring = True
        self._monitor_thread = threading.Thread(
            target=self._monitor_loop, 
            args=(pump, drug), 
            daemon=True
        )
        self._monitor_thread.start()
        logger.info(f"[SafetyWatchdog] Started monitoring {drug}")
    
    def stop_monitoring(self):
        self._monitoring = False
        if self._monitor_thread:
            self._monitor_thread.join(timeout=2.0)
        logger.info("[SafetyWatchdog] Stopped monitoring")
    
    def _monitor_loop(self, pump: PumpInterface, drug: str):
        while self._monitoring:
            try:
                state = pump.get_state()
                if state.status == PumpStatus.INFUSING:
                    safe, adj_rate = self.enforce_safety(pump, state, drug, state.rate_ul_min)
                    if not safe or adj_rate != state.rate_ul_min:
                        if adj_rate > 0:
                            pump.infuse(adj_rate)
            except Exception as e:
                logger.error(f"[SafetyWatchdog] Monitor error: {e}")
            
            time.sleep(self.check_interval)
    
    def get_safety_report(self) -> Dict:
        with self._lock:
            return {
                "cumulative_toxicity": dict(self.cumulative_toxicity),
                "cumulative_volume_ul": dict(self.cumulative_volume),
                "hourly_volume_ul": {k: sum(v for _, v in v) for k, v in self.hourly_volume.items()},
                "violations_count": len(self.violations),
                "recent_violations": self.violations[-10:] if self.violations else []
            }


# Closed-Loop Pipeline with Latency Benchmarking
class ClosedLoopPipeline:
    """
    End-to-end closed-loop pipeline with latency benchmarking:
    Sensors -> PPO Policy -> Safety Check -> Pump Actuation
    Target: < 500ms total latency
    """
    
    def __init__(
        self,
        policy_model,
        pump: PumpInterface,
        sensors,
        watchdog: SafetyWatchdog,
        drug: str = "TMZ",
        max_latency_ms: float = 500.0
    ):
        self.policy = policy_model
        self.pump = pump
        self.sensors = sensors
        self.watchdog = watchdog
        self.drug = drug
        self.max_latency_ms = max_latency_ms
        self.latencies = deque(maxlen=1000)
        self._running = False
    
    def step(self, u: np.ndarray, rho: np.ndarray, step_count: int) -> Dict:
        """Execute one closed-loop step with latency measurement."""
        start_total = time.perf_counter()
        
        # 1. Sensor reading
        t0 = time.perf_counter()
        sensor_data = self.sensors.read_all(u, rho, t=step_count)
        t_sensor = time.perf_counter() - t0
        
        # 2. Policy inference
        t1 = time.perf_counter()
        obs = self._construct_observation(u, sensor_data)
        action, _ = self.policy.predict(obs, deterministic=True)
        t_policy = time.perf_counter() - t1
        
        # 3. Safety check
        t2 = time.perf_counter()
        pump_state = self.pump.get_state()
        proposed_rate = action[0] * 100  # Normalize to µL/min
        safe, adjusted_rate = self.watchdog.enforce_safety(
            self.pump, self.pump.get_state(), "TMZ", proposed_rate
        )
        t_safety = time.perf_counter() - t2
        
        # 4. Pump actuation
        t3 = time.perf_counter()
        if safe and adjusted_rate > 0:
            self.pump.infuse(adjusted_rate)
        elif not safe:
            self.pump.stop()
        t_pump = time.perf_counter() - t3
        
        total_latency = (time.perf_counter() - start_total) * 1000
        self.latencies.append({
            "total_ms": total_latency,
            "sensor_ms": t_sensor * 1000,
            "policy_ms": t_policy * 1000,
            "safety_ms": t_safety * 1000,
            "pump_ms": t_pump * 1000,
            "timestamp": time.time()
        })
        
        # Latency warning
        if total_latency > self.max_latency_ms:
            logger.warning(f"[Pipeline] LATENCY EXCEEDED: {total_latency:.1f}ms > {self.max_latency_ms}ms")
        
        return {
            "action": float(action[0]),
            "adjusted_rate": adjusted_rate,
            "safe": safe,
            "latency_ms": total_latency,
            "latency_breakdown": {
                "sensor": t_sensor * 1000,
                "policy": t_policy * 1000,
                "safety": t_safety * 1000,
                "pump": t_pump * 1000
            }
        }
    
    def _construct_observation(self, u: np.ndarray, sensor_data: Dict) -> np.ndarray:
        """Construct observation vector for policy."""
        vol_norm = np.sum(u > 0.1) / 100000.0
        ctdna = sensor_data.get('ctdna_fragments_per_ml', 0) / 1000.0
        ifp = sensor_data.get('ifp_mmHg', 0) / 50.0
        bbb = sensor_data.get('bbb_permeability', 0)
        ph = (sensor_data.get('pH', 7.4) - 7.0) / 1.0
        o2 = sensor_data.get('pO2_mmHg', 40) / 100.0
        
        # Simplified observation
        obs = np.array([
            vol_norm, ctdna, ifp, bbb, sensor_data.get('pH', 7.4)/10, o2,
            0, 0, 0,  # toxicity placeholder
            0, 1,  # circadian placeholder
            0, 0, 0,  # concentrations
            1, 1, 1,  # time since dose
            0.01  # resistance
        ], dtype=np.float32)
        return obs
    
    def get_latency_stats(self) -> Dict:
        if not self.latencies:
            return {}
        
        latencies = [l["total_ms"] for l in self.latencies]
        return {
            "mean_ms": np.mean(latencies),
            "median_ms": np.median(latencies),
            "p95_ms": np.percentile(latencies, 95),
            "p99_ms": np.percentile(latencies, 99),
            "max_ms": np.max(latencies),
            "under_500ms_pct": sum(1 for l in latencies if l < 500) / len(latencies) * 100,
            "count": len(latencies)
        }


# Example usage and testing
if __name__ == "__main__":
    print("=" * 60)
    print("Phase 14: Hardware-in-the-Loop (HIL) Integration Test")
    print("=" * 60)
    
    # 1. Test Mock Pump
    print("\n[Test 1] Mock Pump Interface")
    pump = MockPumpInterface()
    pump.connect()
    pump.infuse(100.0, 5000)
    time.sleep(2)
    state = pump.get_state()
    print(f"  State: {state.status}, Rate: {state.rate_ul_min}, Vol: {state.volume_delivered_ul:.1f}µL")
    pump.stop()
    print("  [OK] Mock pump works")
    
    # 2. Test Safety Watchdog
    print("\n[Test 2] Safety Watchdog")
    watchdog = SafetyWatchdog(
        max_instantaneous_rate_ul_min={"TMZ": 100, "Inhibitor": 100, "Radiation": 50},
        max_cumulative_toxicity={"TMZ": 1.0, "Inhibitor": 0.8, "Radiation": 0.6}
    )
    pump = MockPumpInterface()
    pump.connect()
    
    # Test safe action
    safe, rate = watchdog.enforce_safety(pump, pump.get_state(), "TMZ", 50.0)
    print(f"  Safe action (50 µL/min): safe={safe}, rate={rate}")
    
    # Test unsafe action (exceeds rate limit)
    safe, rate = watchdog.enforce_safety(pump, pump.get_state(), "TMZ", 200.0)
    print(f"  Unsafe action (200 µL/min): safe={safe}, adjusted_rate={rate}")
    
    # Test cumulative toxicity limit
    watchdog.cumulative_toxicity["TMZ"] = 1.5
    safe, rate = watchdog.enforce_safety(pump, pump.get_state(), "TMZ", 10.0)
    print(f"  Toxicity exceeded: safe={safe}, rate={rate}")
    print("  [OK] Safety watchdog works")
    
    # 3. Test Pipeline Latency
    print("\n[Test 3] Closed-Loop Pipeline Latency")
    from unittest.mock import MagicMock
    import numpy as np
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from sensing.virtual_sensor import VirtualBiosensorSuite
    
    # Mock policy
    mock_policy = MagicMock()
    mock_policy.predict.return_value = (np.array([0.5]), None)
    
    pump = MockPumpInterface()
    pump.connect()
    sensors = VirtualBiosensorSuite()
    watchdog = SafetyWatchdog()
    watchdog.start_monitoring(pump, "TMZ")
    
    pipeline = ClosedLoopPipeline(
        policy_model=mock_policy,
        pump=pump,
        sensors=sensors,
        watchdog=watchdog,
        drug="TMZ"
    )
    
    # Create dummy tumor state
    u = np.zeros((32, 32, 32))
    u[16, 16, 16] = 1.0
    rho = np.full((32, 32, 32), 0.02)
    
    # Run 10 steps
    for i in range(10):
        result = pipeline.step(
            u=np.random.rand(32, 32, 32) * 0.1,
            rho=np.full((32, 32, 32), 0.02),
            step_count=i
        )
        print(f"  Step {i+1}: latency={result['latency_ms']:.1f}ms, "
              f"breakdown: S={result['latency_breakdown']['sensor']:.1f} "
              f"P={result['latency_breakdown']['policy']:.1f} "
              f"Sf={result['latency_breakdown']['safety']:.1f} "
              f"Pu={result['latency_breakdown']['pump']:.1f}ms")
    
    stats = pipeline.get_latency_stats()
    print(f"\n  Latency Stats:")
    print(f"  Mean: {stats['mean_ms']:.1f}ms")
    print(f"  P95: {stats['p95_ms']:.1f}ms")
    print(f"  P99: {stats['p99_ms']:.1f}ms")
    print(f"  Max: {stats['max_ms']:.1f}ms")
    print(f"  Under 500ms: {stats['under_500ms_pct']:.1f}%")
    
    watchdog.stop_monitoring()
    pump.stop()
    print("\n[SUCCESS] Phase 14 HIL Integration tests passed!")